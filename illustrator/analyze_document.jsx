/**********************************************************************
 * analyze_document.jsx
 *
 * Deterministic, READ-ONLY Adobe Illustrator document analyzer.
 *
 * Part of a local scientific-figure normalization system. This script
 * inspects app.activeDocument and exports a structured JSON manifest
 * (plus a human-readable text summary) describing every artboard,
 * layer, group and page item, with rich detail for TextFrameItem,
 * PathItem, RasterItem and PlacedItem objects.
 *
 * DESIGN CONSTRAINTS (why the code looks the way it does):
 *   - Illustrator ExtendScript runs an old (ES3-era) JavaScript engine.
 *     There is NO global JSON object, and modern Array/String helpers
 *     (forEach, map, indexOf, trim, ...) are unreliable. We therefore
 *     hand-roll JSON encoding and use only classic for-loops.
 *   - Many DOM properties are version-dependent or simply absent on
 *     some item subclasses (e.g. PageItem.uuid exists only since
 *     Illustrator 24.0 / 2020; controlBounds may be missing; there is
 *     NO generic rotation property on page items). Every property read
 *     goes through a safe wrapper that records a warning and continues.
 *   - The document is never modified. We only read.
 *
 * USAGE:  File > Scripts > Other Script...  (see README.md)
 *
 * Author: scientific-figure normalization system
 *********************************************************************/

#target illustrator

(function () {
    "use strict";

    /* ============================================================== *
     *  0. Constants / metadata
     * ============================================================== */

    var SCRIPT_NAME    = "analyze_document.jsx";
    var SCHEMA_VERSION = "1.0.0";

    /* ============================================================== *
     *  1. Shared collector state
     * ============================================================== */

    // Global accumulators. Kept in one object so helper functions can
    // reach them without a web of parameters.
    var C = {
        warnings:   [],   // array of { context, message }
        typeCounts: {},   // typename -> integer count
        fontUsage:  {},   // "Family / Style" -> { name, family, style, characters, runs }
        totalItems: 0,
        maxDepth:   0
    };

    /* --- warning + type counting helpers --------------------------- */

    function addWarning(context, message) {
        C.warnings.push({ context: String(context), message: String(message) });
    }

    function countType(typename) {
        var key = String(typename);
        if (C.typeCounts[key] === undefined) { C.typeCounts[key] = 0; }
        C.typeCounts[key] += 1;
        C.totalItems += 1;
    }

    function noteFont(name, family, style, nChars) {
        var key = String(family) + " / " + String(style);
        var rec = C.fontUsage[key];
        if (!rec) {
            rec = { name: String(name), family: String(family),
                    style: String(style), characters: 0, runs: 0 };
            C.fontUsage[key] = rec;
        }
        rec.runs += 1;
        rec.characters += nChars;
    }

    /* ============================================================== *
     *  2. Safe property access
     *
     *  Every DOM read that might throw or be unsupported is funneled
     *  through here so a single bad property never aborts the analysis.
     * ============================================================== */

    // Run fn(); on any exception, log a warning and return dflt.
    function safe(fn, context, dflt) {
        try {
            var v = fn();
            if (v === undefined) { return dflt; }
            return v;
        } catch (e) {
            addWarning(context, e && e.message ? e.message : String(e));
            return dflt;
        }
    }

    // Convert a value to a finite number or null.
    function num(v) {
        if (v === null || v === undefined) { return null; }
        var n = Number(v);
        if (isNaN(n) || !isFinite(n)) { return null; }
        return n;
    }

    // Round to a fixed number of decimals to reduce float noise while
    // keeping geometry meaningful.
    function round(v, dp) {
        var n = num(v);
        if (n === null) { return null; }
        var f = Math.pow(10, dp === undefined ? 4 : dp);
        return Math.round(n * f) / f;
    }

    // Copy an Illustrator Rect ([left, top, right, bottom]) into a
    // plain, rounded 4-element array. Returns null if unavailable.
    function boundsToArray(rectFn, context) {
        return safe(function () {
            var r = rectFn();
            if (!r || r.length < 4) { return null; }
            return [round(r[0]), round(r[1]), round(r[2]), round(r[3])];
        }, context, null);
    }

    // Stringify an ExtendScript enum (e.g. TextType.POINTTEXT) safely.
    function enumStr(fn, context) {
        return safe(function () {
            var v = fn();
            if (v === null || v === undefined) { return null; }
            return String(v);
        }, context, null);
    }

    /* ============================================================== *
     *  3. Hand-rolled JSON encoder (ES3-safe, pretty printed)
     * ============================================================== */

    function pad4(hex) {
        while (hex.length < 4) { hex = "0" + hex; }
        return hex;
    }

    // Escape a string into a valid JSON string literal WITHOUT relying
    // on JSON.stringify (unavailable in ExtendScript).
    function jsonString(s) {
        s = String(s);
        var out = '"';
        for (var i = 0; i < s.length; i++) {
            var c    = s.charAt(i);
            var code = s.charCodeAt(i);
            if (c === '"')       { out += '\\"'; }
            else if (c === '\\') { out += '\\\\'; }
            else if (c === '\b') { out += '\\b'; }
            else if (c === '\f') { out += '\\f'; }
            else if (c === '\n') { out += '\\n'; }
            else if (c === '\r') { out += '\\r'; }
            else if (c === '\t') { out += '\\t'; }
            else if (code < 0x20 || code === 0x2028 || code === 0x2029) {
                out += '\\u' + pad4(code.toString(16));
            } else {
                out += c;
            }
        }
        return out + '"';
    }

    // Recursively encode a plain value (null/number/boolean/string/
    // Array/Object) into pretty-printed JSON. Object key order follows
    // insertion order, which ExtendScript preserves.
    function jsonEncode(value, indent) {
        if (indent === undefined) { indent = ""; }
        var childIndent = indent + "  ";

        if (value === null || value === undefined) { return "null"; }

        var t = typeof value;
        if (t === "number")  { return isFinite(value) ? String(value) : "null"; }
        if (t === "boolean") { return value ? "true" : "false"; }
        if (t === "string")  { return jsonString(value); }

        if (value instanceof Array) {
            if (value.length === 0) { return "[]"; }
            var items = [];
            for (var i = 0; i < value.length; i++) {
                items.push(childIndent + jsonEncode(value[i], childIndent));
            }
            return "[\n" + items.join(",\n") + "\n" + indent + "]";
        }

        // Plain object
        var keys = [];
        for (var k in value) {
            if (value.hasOwnProperty(k)) { keys.push(k); }
        }
        if (keys.length === 0) { return "{}"; }
        var parts = [];
        for (var j = 0; j < keys.length; j++) {
            var key = keys[j];
            parts.push(childIndent + jsonString(key) + ": " +
                       jsonEncode(value[key], childIndent));
        }
        return "{\n" + parts.join(",\n") + "\n" + indent + "}";
    }

    /* ============================================================== *
     *  4. Color serialization
     *
     *  Colors may be RGB / CMYK / Gray / Spot / Lab / Pattern /
     *  Gradient / NoColor. We switch on .typename and read only the
     *  fields valid for that subtype, each defensively.
     * ============================================================== */

    function serializeColor(colorFn, context) {
        var color = safe(colorFn, context, null);
        if (!color) { return null; }

        var typename = safe(function () { return color.typename; },
                            context + ".typename", "UnknownColor");
        var out = { type: String(typename) };

        switch (String(typename)) {
        case "RGBColor":
            out.red   = round(safe(function () { return color.red;   }, context, null), 3);
            out.green = round(safe(function () { return color.green; }, context, null), 3);
            out.blue  = round(safe(function () { return color.blue;  }, context, null), 3);
            break;
        case "CMYKColor":
            out.cyan    = round(safe(function () { return color.cyan;    }, context, null), 3);
            out.magenta = round(safe(function () { return color.magenta; }, context, null), 3);
            out.yellow  = round(safe(function () { return color.yellow;  }, context, null), 3);
            out.black   = round(safe(function () { return color.black;   }, context, null), 3);
            break;
        case "GrayColor":
            out.gray = round(safe(function () { return color.gray; }, context, null), 3);
            break;
        case "LabColor":
            out.l = round(safe(function () { return color.l; }, context, null), 3);
            out.a = round(safe(function () { return color.a; }, context, null), 3);
            out.b = round(safe(function () { return color.b; }, context, null), 3);
            break;
        case "SpotColor":
            out.spotName = safe(function () { return color.spot.name; }, context, null);
            out.tint     = round(safe(function () { return color.tint; }, context, null), 3);
            break;
        case "GradientColor":
            out.gradientName = safe(function () { return color.gradient.name; }, context, null);
            break;
        case "PatternColor":
            out.patternName = safe(function () { return color.pattern.name; }, context, null);
            break;
        case "NoColor":
            // Nothing further to record.
            break;
        default:
            addWarning(context, "Unhandled color typename: " + typename);
            break;
        }
        return out;
    }

    /* ============================================================== *
     *  5. Generic page-item properties (common to every item)
     * ============================================================== */

    function commonItemInfo(item, parentDesc) {
        var typename = safe(function () { return item.typename; },
                            "item.typename", "PageItem");
        countType(typename);

        var info = {};
        info.typename = String(typename);

        // uuid: only since Illustrator 24.0 (2020). May throw on older.
        info.uuid = safe(function () { return item.uuid; }, "item.uuid", null);

        info.name = safe(function () { return item.name; }, "item.name", "");
        info.note = safe(function () { return item.note; }, "item.note", "");

        // State flags.
        info.locked   = safe(function () { return item.locked;   }, "item.locked",   null);
        info.hidden   = safe(function () { return item.hidden;   }, "item.hidden",   null);
        // Not all subclasses expose .editable; read defensively.
        info.editable = safe(function () { return item.editable; }, "item.editable", null);

        // Bounds. controlBounds may be absent on some versions/items.
        info.geometricBounds = boundsToArray(function () { return item.geometricBounds; }, "item.geometricBounds");
        info.visibleBounds   = boundsToArray(function () { return item.visibleBounds;   }, "item.visibleBounds");
        info.controlBounds   = boundsToArray(function () { return item.controlBounds;   }, "item.controlBounds");

        // Position / size.
        info.position = safe(function () {
            var p = item.position;
            if (!p || p.length < 2) { return null; }
            return [round(p[0]), round(p[1])];
        }, "item.position", null);
        info.width  = round(safe(function () { return item.width;  }, "item.width",  null));
        info.height = round(safe(function () { return item.height; }, "item.height", null));

        // Rotation: the Illustrator DOM does NOT expose a generic
        // rotation angle on page items (rotation lives inside the
        // transformation matrix, which is not readable). Reported as
        // null; a single global warning is emitted in main().
        info.rotation = null;

        // Parent hierarchy description (shallow, to avoid cycles).
        info.parent = parentDesc ? parentDesc : null;

        return info;
    }

    // Build a short, non-recursive descriptor of a container for the
    // child's "parent" field.
    function describeParent(container, kind) {
        return {
            kind:     kind,
            typename: safe(function () { return container.typename; }, "parent.typename", kind),
            name:     safe(function () { return container.name; }, "parent.name", "")
        };
    }

    /* ============================================================== *
     *  6. TextFrameItem detail
     * ============================================================== */

    // Build a signature string for a character's formatting so that
    // consecutive characters with identical formatting collapse into a
    // single run.
    function runSignature(fontName, size, baseline, hScale, vScale,
                          tracking, caps, fill, stroke) {
        return [fontName, size, baseline, hScale, vScale, tracking, caps,
                jsonEncode(fill), jsonEncode(stroke)].join("|");
    }

    // Read one character's formatting into a plain object.
    function readCharFormatting(ca) {
        var fontName   = safe(function () { return ca.textFont.name;   }, "char.textFont.name",   null);
        var fontFamily = safe(function () { return ca.textFont.family; }, "char.textFont.family", null);
        var fontStyle  = safe(function () { return ca.textFont.style;  }, "char.textFont.style",  null);

        return {
            fontName:        fontName,
            fontFamily:      fontFamily,
            fontStyle:       fontStyle,
            size:            round(safe(function () { return ca.size;            }, "char.size",            null), 3),
            baselineShift:   round(safe(function () { return ca.baselineShift;   }, "char.baselineShift",   null), 3),
            horizontalScale: round(safe(function () { return ca.horizontalScale; }, "char.horizontalScale", null), 3),
            verticalScale:   round(safe(function () { return ca.verticalScale;   }, "char.verticalScale",   null), 3),
            tracking:        round(safe(function () { return ca.tracking;        }, "char.tracking",        null), 3),
            capitalization:  enumStr(function () { return ca.capitalization; }, "char.capitalization"),
            fillColor:       serializeColor(function () { return ca.fillColor;   }, "char.fillColor"),
            strokeColor:     serializeColor(function () { return ca.strokeColor; }, "char.strokeColor")
        };
    }

    // Walk the characters of a text frame and coalesce them into runs
    // of identical formatting.
    function extractRuns(tf) {
        var runs  = [];
        var chars = safe(function () { return tf.textRange.characters; }, "tf.textRange.characters", null);
        if (!chars) { return runs; }

        var len = safe(function () { return chars.length; }, "tf.characters.length", 0);

        var curr = null;      // current run object being built
        var currSig = null;   // its signature

        for (var i = 0; i < len; i++) {
            var ch = safe(function () { return chars[i]; }, "char[" + i + "]", null);
            if (!ch) { continue; }

            var content = safe(function () { return ch.contents; }, "char[" + i + "].contents", "");
            var ca = safe(function () { return ch.characterAttributes; }, "char[" + i + "].attrs", null);
            if (!ca) { continue; }

            var fmt = readCharFormatting(ca);
            var sig = runSignature(fmt.fontName, fmt.size, fmt.baselineShift,
                                   fmt.horizontalScale, fmt.verticalScale,
                                   fmt.tracking, fmt.capitalization,
                                   fmt.fillColor, fmt.strokeColor);

            if (curr && sig === currSig) {
                curr.text += content;
                curr.length += 1;
            } else {
                if (curr) { runs.push(curr); }
                curr = {
                    text:            content,
                    length:          1,
                    font:            { name: fmt.fontName, family: fmt.fontFamily, style: fmt.fontStyle },
                    size:            fmt.size,
                    baselineShift:   fmt.baselineShift,
                    horizontalScale: fmt.horizontalScale,
                    verticalScale:   fmt.verticalScale,
                    tracking:        fmt.tracking,
                    capitalization:  fmt.capitalization,
                    fillColor:       fmt.fillColor,
                    strokeColor:     fmt.strokeColor
                };
                currSig = sig;
            }
        }
        if (curr) { runs.push(curr); }

        // Record font usage statistics from the finished runs.
        for (var r = 0; r < runs.length; r++) {
            var f = runs[r].font;
            if (f && f.name) {
                noteFont(f.name, f.family, f.style, runs[r].length);
            }
        }
        return runs;
    }

    // Per-paragraph justification list.
    function extractParagraphs(tf) {
        var out = [];
        var paras = safe(function () { return tf.paragraphs; }, "tf.paragraphs", null);
        if (!paras) { return out; }
        var n = safe(function () { return paras.length; }, "tf.paragraphs.length", 0);
        for (var i = 0; i < n; i++) {
            out.push({
                index: i,
                justification: enumStr(function () {
                    return paras[i].paragraphAttributes.justification;
                }, "paragraph[" + i + "].justification")
            });
        }
        return out;
    }

    function textFrameDetail(tf) {
        var d = {};
        d.contents    = safe(function () { return tf.contents; },    "tf.contents", "");
        d.kind        = enumStr(function () { return tf.kind; },        "tf.kind");
        d.orientation = enumStr(function () { return tf.orientation; }, "tf.orientation");
        d.characterCount = safe(function () { return tf.textRange.characters.length; },
                                "tf.characterCount", null);
        d.paragraphs = extractParagraphs(tf);
        d.runs       = extractRuns(tf);
        return d;
    }

    /* ============================================================== *
     *  7. PathItem detail
     * ============================================================== */

    function pathItemDetail(p) {
        var d = {};
        d.closed        = safe(function () { return p.closed;  }, "path.closed",  null);
        d.filled        = safe(function () { return p.filled;  }, "path.filled",  null);
        d.stroked       = safe(function () { return p.stroked; }, "path.stroked", null);
        d.clipping      = safe(function () { return p.clipping; }, "path.clipping", null);
        d.strokeWidth   = round(safe(function () { return p.strokeWidth; }, "path.strokeWidth", null), 4);
        d.pathPointCount = safe(function () { return p.pathPoints.length; }, "path.pathPoints.length", null);
        d.fillColor     = serializeColor(function () { return p.fillColor;   }, "path.fillColor");
        d.strokeColor   = serializeColor(function () { return p.strokeColor; }, "path.strokeColor");
        return d;
    }

    /* ============================================================== *
     *  8. Raster / Placed detail
     * ============================================================== */

    function fileInfo(fileFn, context) {
        var f = safe(fileFn, context, null);
        if (!f) { return null; }
        return {
            fsName: safe(function () { return f.fsName; }, context + ".fsName", null),
            exists: safe(function () { return f.exists; }, context + ".exists", null)
        };
    }

    function rasterItemDetail(ri) {
        var d = {};
        // RasterItem.embedded indicates whether pixels are embedded.
        d.embedded = safe(function () { return ri.embedded; }, "raster.embedded", null);
        d.file     = fileInfo(function () { return ri.file; }, "raster.file");
        d.linked   = (d.embedded === false);
        d.effectiveBounds = boundsToArray(function () { return ri.boundingBox; }, "raster.boundingBox");
        if (d.effectiveBounds === null) {
            d.effectiveBounds = boundsToArray(function () { return ri.geometricBounds; }, "raster.geometricBounds");
        }
        return d;
    }

    function placedItemDetail(pi) {
        var d = {};
        // A PlacedItem is a linked/placed artwork reference. It carries
        // a .file to the external asset; when embedded it is converted
        // to a RasterItem, so we report embedded=false here.
        d.embedded = false;
        d.linked   = true;
        d.file     = fileInfo(function () { return pi.file; }, "placed.file");
        d.effectiveBounds = boundsToArray(function () { return pi.geometricBounds; }, "placed.geometricBounds");
        return d;
    }

    /* ============================================================== *
     *  9. Recursive item walker
     *
     *  Dispatches on typename, fills in type-specific detail, and
     *  recurses into containers (GroupItem, CompoundPathItem).
     * ============================================================== */

    function processItem(item, parentDesc, depth) {
        if (depth > C.maxDepth) { C.maxDepth = depth; }

        var info = commonItemInfo(item, parentDesc);
        var typename = info.typename;

        // IMPORTANT — runtime typename vs. scripting class name:
        // A text frame's scripting CLASS is documented as "TextFrameItem",
        // but the value returned by PageItem.typename at runtime is the
        // string "TextFrame" (Adobe's own samples test `typename == "TextFrame"`).
        // We accept both: "TextFrame" is what real documents actually report
        // and is required; "TextFrameItem" is kept as a defensive alias.
        //
        // Class-name/typename review (verified against a 1108-object
        // production document): PathItem, CompoundPathItem, GroupItem and
        // RasterItem all report typename == class name and dispatch correctly
        // below. PlacedItem likewise reports "PlacedItem". "TextFrame" is the
        // only page-item class whose runtime typename differs from its class
        // name, so it is the only case needing this alias.
        if (typename === "TextFrame" || typename === "TextFrameItem") {
            info.text = textFrameDetail(item);

        } else if (typename === "PathItem") {
            info.path = pathItemDetail(item);

        } else if (typename === "RasterItem") {
            info.raster = rasterItemDetail(item);

        } else if (typename === "PlacedItem") {
            info.placed = placedItemDetail(item);

        } else if (typename === "CompoundPathItem") {
            // A compound path contains child PathItems (holes, etc.).
            info.compound = { closed: safe(function () { return item.pathItems.length; },
                                            "compound.pathItems.length", null) };
            info.children = processChildren(item, describeParent(item, "compoundPath"), depth + 1);

        } else if (typename === "GroupItem") {
            info.group = {
                clipped: safe(function () { return item.clipped; }, "group.clipped", null)
            };
            info.children = processChildren(item, describeParent(item, "group"), depth + 1);

        } else {
            // SymbolItem, MeshItem, GraphItem, PluginItem, etc.:
            // captured with common info only. Record a note so callers
            // know detail was not specialized.
            addWarning("processItem", "No specialized handler for typename '" +
                       typename + "' (captured with generic fields only).");
        }

        return info;
    }

    // Iterate the DIRECT children of a container via .pageItems.
    // (In the Illustrator DOM, container.pageItems yields immediate
    // children; nested groups are recursed explicitly above.)
    function processChildren(container, parentDesc, depth) {
        var out = [];
        var items = safe(function () { return container.pageItems; }, "container.pageItems", null);
        if (!items) { return out; }
        var n = safe(function () { return items.length; }, "container.pageItems.length", 0);
        for (var i = 0; i < n; i++) {
            var child = safe(function () { return items[i]; }, "pageItems[" + i + "]", null);
            if (!child) { continue; }
            out.push(processItem(child, parentDesc, depth));
        }
        return out;
    }

    /* ============================================================== *
     * 10. Layer walker (handles nested sublayers)
     * ============================================================== */

    function processLayer(layer, depth) {
        var node = { kind: "layer" };
        node.name     = safe(function () { return layer.name; },     "layer.name", "");
        node.visible  = safe(function () { return layer.visible; },  "layer.visible", null);
        node.locked   = safe(function () { return layer.locked; },   "layer.locked", null);
        node.opacity  = round(safe(function () { return layer.opacity; }, "layer.opacity", null), 2);
        node.color    = serializeColor(function () { return layer.color; }, "layer.color");

        var parentDesc = describeParent(layer, "layer");

        // Direct page items of this layer.
        node.items = processChildren(layer, parentDesc, depth + 1);

        // Nested sublayers.
        node.sublayers = [];
        var subs = safe(function () { return layer.layers; }, "layer.layers", null);
        if (subs) {
            var n = safe(function () { return subs.length; }, "layer.layers.length", 0);
            for (var i = 0; i < n; i++) {
                var sub = safe(function () { return subs[i]; }, "layer.layers[" + i + "]", null);
                if (sub) { node.sublayers.push(processLayer(sub, depth + 1)); }
            }
        }
        return node;
    }

    /* ============================================================== *
     * 11. Document / artboard metadata
     * ============================================================== */

    function documentInfo(doc) {
        var d = {};
        d.name     = safe(function () { return doc.name; }, "doc.name", "");
        d.fullPath = safe(function () {
            return doc.fullName ? doc.fullName.fsName : null;
        }, "doc.fullName", null);
        d.saved     = safe(function () { return doc.saved; }, "doc.saved", null);
        d.colorSpace = enumStr(function () { return doc.documentColorSpace; }, "doc.documentColorSpace");
        d.rulerUnits = enumStr(function () { return doc.rulerUnits; }, "doc.rulerUnits");

        // Illustrator application version.
        d.illustratorVersion = safe(function () { return app.version; }, "app.version", null);
        d.buildNumber        = safe(function () { return app.buildNumber; }, "app.buildNumber", null);

        // Artboards.
        d.artboards = [];
        var abs = safe(function () { return doc.artboards; }, "doc.artboards", null);
        if (abs) {
            var n = safe(function () { return abs.length; }, "doc.artboards.length", 0);
            for (var i = 0; i < n; i++) {
                (function (idx) {
                    var ab = safe(function () { return abs[idx]; }, "artboard[" + idx + "]", null);
                    if (!ab) { return; }
                    var rect = boundsToArray(function () { return ab.artboardRect; }, "artboard.rect");
                    var w = null, h = null;
                    if (rect) {
                        w = round(rect[2] - rect[0]);       // right - left
                        h = round(rect[1] - rect[3]);       // top - bottom
                    }
                    d.artboards.push({
                        index: idx,
                        name: safe(function () { return ab.name; }, "artboard.name", ""),
                        artboardRect: rect,
                        width: w,
                        height: h
                    });
                })(i);
            }
        }
        return d;
    }

    /* ============================================================== *
     * 12. Statistics / summary assembly
     * ============================================================== */

    function buildStatistics() {
        // Sorted type-count list for stable output.
        var typeKeys = [];
        for (var k in C.typeCounts) {
            if (C.typeCounts.hasOwnProperty(k)) { typeKeys.push(k); }
        }
        typeKeys.sort();
        var byType = [];
        for (var i = 0; i < typeKeys.length; i++) {
            byType.push({ typename: typeKeys[i], count: C.typeCounts[typeKeys[i]] });
        }

        // Sorted font-usage list.
        var fontKeys = [];
        for (var fk in C.fontUsage) {
            if (C.fontUsage.hasOwnProperty(fk)) { fontKeys.push(fk); }
        }
        fontKeys.sort();
        var fonts = [];
        for (var j = 0; j < fontKeys.length; j++) {
            fonts.push(C.fontUsage[fontKeys[j]]);
        }

        return {
            totalItems: C.totalItems,
            maxNestingDepth: C.maxDepth,
            itemsByType: byType,
            fonts: fonts,
            warningCount: C.warnings.length
        };
    }

    // Human-readable text summary (CRLF for Windows friendliness).
    function buildTextSummary(manifest) {
        var NL = "\r\n";
        var s = [];
        s.push("Illustrator Document Analysis Summary");
        s.push("=====================================");
        s.push("Generated by : " + SCRIPT_NAME + " (schema " + SCHEMA_VERSION + ")");
        s.push("");
        s.push("Document     : " + manifest.document.name);
        s.push("Full path    : " + (manifest.document.fullPath || "(unsaved / no path)"));
        s.push("AI version   : " + manifest.document.illustratorVersion);
        s.push("Color space  : " + manifest.document.colorSpace);
        s.push("Artboards    : " + manifest.document.artboards.length);
        s.push("");
        s.push("Object counts");
        s.push("-------------");
        var t = manifest.statistics.itemsByType;
        for (var i = 0; i < t.length; i++) {
            s.push("  " + rpad(t[i].typename, 20) + t[i].count);
        }
        s.push("  " + rpad("TOTAL", 20) + manifest.statistics.totalItems);
        s.push("  max nesting depth : " + manifest.statistics.maxNestingDepth);
        s.push("");
        s.push("Font usage (family / style : characters in N runs)");
        s.push("--------------------------------------------------");
        var f = manifest.statistics.fonts;
        if (f.length === 0) {
            s.push("  (no text found)");
        } else {
            for (var k = 0; k < f.length; k++) {
                s.push("  " + rpad(f[k].family + " / " + f[k].style, 32) +
                       f[k].characters + " chars in " + f[k].runs + " runs" +
                       "   [" + f[k].name + "]");
            }
        }
        s.push("");
        s.push("Warnings     : " + manifest.warnings.length);
        if (manifest.warnings.length > 0) {
            s.push("  (see JSON manifest 'warnings' array for full detail)");
            var shown = manifest.warnings.length < 20 ? manifest.warnings.length : 20;
            for (var w = 0; w < shown; w++) {
                s.push("  - [" + manifest.warnings[w].context + "] " +
                       manifest.warnings[w].message);
            }
            if (manifest.warnings.length > shown) {
                s.push("  ... and " + (manifest.warnings.length - shown) + " more.");
            }
        }
        s.push("");
        return s.join(NL);
    }

    function rpad(str, width) {
        str = String(str);
        while (str.length < width) { str += " "; }
        return str;
    }

    /* ============================================================== *
     * 13. File output (UTF-8, overwrite-protected)
     * ============================================================== */

    // Write text to a File as UTF-8. Returns true on success.
    function writeUTF8(file, text) {
        file.encoding = "UTF-8";
        var opened = file.open("w");
        if (!opened) {
            addWarning("writeUTF8", "Could not open for writing: " + file.fsName);
            return false;
        }
        file.write(text);
        file.close();
        return true;
    }

    // Refuse to overwrite an existing file unless the user confirms.
    function confirmWrite(file) {
        if (file.exists) {
            return confirm("File already exists:\n" + file.fsName +
                           "\n\nOverwrite it?");
        }
        return true;
    }

    /* ============================================================== *
     * 14. Main
     * ============================================================== */

    function main() {
        // Guard: an active document is required.
        if (app.documents.length === 0) {
            alert("No document is open. Open an Illustrator document and run again.");
            return;
        }
        var doc = app.activeDocument;

        // One-time note about the rotation limitation.
        addWarning("rotation",
            "The Illustrator DOM does not expose a per-item rotation angle; " +
            "'rotation' is reported as null for all page items.");

        // Ask the user where to write output. Read-only up to this point.
        var outDir = Folder.selectDialog("Choose an output directory for the analysis manifest");
        if (!outDir) {
            // User cancelled: do nothing, modify nothing.
            return;
        }

        // Assemble the manifest tree.
        var manifest = {};
        manifest.schemaVersion = SCHEMA_VERSION;
        manifest.generator     = SCRIPT_NAME;
        manifest.generatedAt   = (new Date()).toString();
        manifest.document      = documentInfo(doc);

        manifest.layers = [];
        var layers = safe(function () { return doc.layers; }, "doc.layers", null);
        if (layers) {
            var n = safe(function () { return layers.length; }, "doc.layers.length", 0);
            for (var i = 0; i < n; i++) {
                var lyr = safe(function () { return layers[i]; }, "doc.layers[" + i + "]", null);
                if (lyr) { manifest.layers.push(processLayer(lyr, 0)); }
            }
        }

        manifest.statistics = buildStatistics();
        manifest.warnings   = C.warnings;

        // Derive output base name from the document name.
        var docName = manifest.document.name || "untitled";
        var base = docName.replace(/\.[Aa][Ii]$/, "").replace(/\.[^.]+$/, "");
        var jsonFile = new File(outDir.fsName + "/" + base + "_manifest.json");
        var txtFile  = new File(outDir.fsName + "/" + base + "_summary.txt");

        // Overwrite protection for BOTH outputs.
        if (!confirmWrite(jsonFile)) { alert("Aborted: JSON output not overwritten."); return; }
        if (!confirmWrite(txtFile))  { alert("Aborted: text summary not overwritten."); return; }

        var jsonText = jsonEncode(manifest, "");
        var txtText  = buildTextSummary(manifest);

        var okJson = writeUTF8(jsonFile, jsonText);
        var okTxt  = writeUTF8(txtFile, txtText);

        var msg = "Analysis complete.\n\n" +
                  "Items analyzed : " + manifest.statistics.totalItems + "\n" +
                  "Warnings       : " + manifest.warnings.length + "\n\n" +
                  (okJson ? ("JSON   : " + jsonFile.fsName + "\n") : "JSON   : FAILED\n") +
                  (okTxt  ? ("Summary: " + txtFile.fsName)          : "Summary: FAILED");
        alert(msg);
    }

    // Top-level guard so any unexpected error is surfaced, not swallowed.
    try {
        main();
    } catch (e) {
        alert("analyze_document.jsx failed:\n" + (e && e.message ? e.message : String(e)) +
              (e && e.line ? ("\n(line " + e.line + ")") : ""));
    }

    /* ============================================================== *
     * 15. Test-only seam (inert inside Illustrator)
     *
     *  ExtendScript has no CommonJS `module`, so `typeof module` is
     *  "undefined" inside Illustrator and this block never runs there:
     *  it does NOT change analyzer behavior. A Node-based regression
     *  harness supplies a `module` object to reach the pure dispatch and
     *  text-frame handler for testing. See tests/textframe_dispatch.test.js.
     * ============================================================== */
    if (typeof module !== "undefined" && module && module.exports) {
        module.exports = {
            processItem:     processItem,
            textFrameDetail: textFrameDetail,
            getWarnings:     function () { return C.warnings; }
        };
    }

})();
