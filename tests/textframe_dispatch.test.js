/**
 * Regression test for the TextFrame typename-dispatch bug.
 *
 * Background:
 *   Illustrator's scripting CLASS is documented as "TextFrameItem", but the
 *   value returned by PageItem.typename at RUNTIME is the string "TextFrame".
 *   An earlier version of analyze_document.jsx dispatched on "TextFrameItem"
 *   only, so all 216 text frames in Figure2.ai fell through to the generic
 *   handler and no text / fonts / character runs were extracted.
 *
 * This test proves that an object whose typename is "TextFrame" now enters
 * the specialized text-frame handler (info.text populated with contents and
 * character-format runs) and produces NO "No specialized handler" warning.
 *
 * Run:   node tests/textframe_dispatch.test.js
 * Exit:  0 = pass, non-zero = fail.
 *
 * The analyzer is loaded through a small Node "test seam" (module.exports)
 * that is inert inside Illustrator. ExtendScript preprocessor directives
 * (#target/#include/#targetengine) are stripped before evaluation because
 * they are not valid plain JavaScript.
 */

"use strict";

var fs     = require("fs");
var path   = require("path");
var vm     = require("vm");
var assert = require("assert");

/* ---- 1. Load analyzer source, stripping ExtendScript directives -------- */

var srcPath = path.join(__dirname, "..", "illustrator", "analyze_document.jsx");
var src = fs.readFileSync(srcPath, "utf8")
    .split(/\r?\n/)
    .map(function (line) { return /^\s*#/.test(line) ? "" : line; })  // drop #target etc.
    .join("\n");

/* ---- 2. Minimal Illustrator-global sandbox ----------------------------- *
 * app.documents.length === 0 makes main() early-exit cleanly (no open doc),
 * so evaluating the file has no side effects beyond running the test seam.  */

var moduleShim = { exports: {} };
var sandbox = {
    app:     { documents: { length: 0 }, name: "Adobe Illustrator", version: "test" },
    alert:   function () {},
    confirm: function () { return true; },
    Folder:  function () {},
    File:    function () {},
    module:  moduleShim,
    console: console
};
vm.createContext(sandbox);
vm.runInContext(src, sandbox, { filename: "analyze_document.jsx" });

var api = moduleShim.exports;
assert.ok(api && typeof api.processItem === "function",
    "analyzer must export processItem via the test seam");

/* ---- 3. Fake page item whose RUNTIME typename is "TextFrame" ------------ */

function makeChar(ch) {
    return {
        contents: ch,
        characterAttributes: {
            textFont:        { name: "ArialMT", family: "Arial", style: "Regular" },
            size:            8,
            baselineShift:   0,
            horizontalScale: 100,
            verticalScale:   100,
            tracking:        0,
            capitalization:  "FontCapsOption.NORMALCAPS",
            fillColor:       { typename: "RGBColor", red: 0, green: 0, blue: 0 },
            strokeColor:     { typename: "NoColor" }
        }
    };
}

var CONTENTS = "Preferred phase";
var chars = CONTENTS.split("").map(makeChar);   // real array: has .length + [i]

var fakeTextFrame = {
    typename: "TextFrame",              // <-- the runtime value, not the class name
    uuid: "TEST-UUID-0001",
    name: "phase-label", note: "",
    locked: false, hidden: false, editable: true,
    geometricBounds: [100, 200, 92, 150],
    visibleBounds:   [100, 200, 92, 150],
    controlBounds:   [100, 200, 92, 150],
    position: [150, 100], width: 50, height: 8,
    contents: CONTENTS,
    kind: "TextType.POINTTEXT",
    orientation: "TextOrientation.HORIZONTAL",
    paragraphs: [ { paragraphAttributes: { justification: "Justification.LEFT" } } ],
    textRange: { characters: chars }
};
var parent = { kind: "layer", typename: "Layer", name: "Figure" };

/* ---- 4. Dispatch ------------------------------------------------------- */

var info = api.processItem(fakeTextFrame, parent, 0);

/* ---- 5. Assertions: specialized handler ran, generic fallback did not --- */

assert.strictEqual(info.typename, "TextFrame", "typename preserved");
assert.ok(info.text, "info.text must be populated (specialized handler ran)");
assert.strictEqual(info.text.contents, CONTENTS, "contents extracted");
assert.ok(Array.isArray(info.text.runs) && info.text.runs.length >= 1,
    "character runs extracted");
assert.strictEqual(info.text.runs[0].font.family, "Arial", "run font family");
assert.strictEqual(info.text.runs[0].font.style, "Regular", "run font style");
assert.strictEqual(info.text.runs[0].size, 8, "run size (pt)");
assert.strictEqual(info.text.kind, "TextType.POINTTEXT", "text-frame kind captured");

/* ---- 6. No "unsupported TextFrame" warning ----------------------------- */

var warns = api.getWarnings();
var badge = warns.filter(function (w) {
    return /No specialized handler for typename 'TextFrame'/.test(w.message);
});
assert.strictEqual(badge.length, 0,
    "no 'No specialized handler for typename TextFrame' warning should be produced");

/* ---- 7. Guard: a genuinely unknown type still warns (dispatch intact) --- */

var beforeUnknown = api.getWarnings().length;
var fakeMesh = { typename: "MeshItem", name: "", note: "" };
api.processItem(fakeMesh, parent, 0);
var afterUnknown = api.getWarnings().length;
assert.ok(afterUnknown > beforeUnknown,
    "unknown typename 'MeshItem' must still record a warning (fallback preserved)");

console.log("PASS: 'TextFrame' routes to the specialized text handler.");
console.log("      contents=" + JSON.stringify(info.text.contents) +
            "  runs=" + info.text.runs.length +
            "  font=" + info.text.runs[0].font.family + "/" + info.text.runs[0].font.style +
            "  size=" + info.text.runs[0].size + "pt");
