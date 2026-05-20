"""
Video Trim/Crop/Speed Tool
- Visual crop: drag rectangle on video
- Trim: set in/out points with buttons or type timestamps
- Speed: adjust playback speed multiplier
- Export: FFmpeg with audio, PowerPoint-friendly MP4 (H.264 + AAC)
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import cv2
from PIL import Image, ImageTk
import subprocess
import shutil
import os
import threading

class VideoTool:
    def __init__(self, root):
        self.root = root
        self.root.title("Video Trim/Crop/Speed Tool")
        self.root.configure(bg="#2b2b2b")

        self.cap = None
        self.total_frames = 0
        self.fps = 30
        self.duration = 0
        self.width = 0
        self.height = 0
        self.filepath = None

        # Display scaling
        self.scale = 1.0

        # Crop state
        self.crop_start = None
        self.crop_rect = None  # (x1, y1, x2, y2) in original video coords
        self.dragging = False
        self.drag_start = None

        # Trim state
        self.trim_in = 0.0
        self.trim_out = None

        # Playback state
        self.playing = False
        self.current_frame = 0

        self._build_ui()

    def _build_ui(self):
        frame_style = {"bg": "#2b2b2b"}
        label_style = {"bg": "#2b2b2b", "fg": "#ffffff"}
        btn_style = {"bg": "#444444", "fg": "#ffffff", "activebackground": "#666666",
                     "activeforeground": "#ffffff", "relief": "flat", "padx": 8, "pady": 4}

        # Top bar
        top = tk.Frame(self.root, **frame_style)
        top.pack(fill="x", padx=5, pady=5)
        tk.Button(top, text="Open Video", command=self.open_file, **btn_style).pack(side="left")
        self.file_label = tk.Label(top, text="No file loaded", **label_style)
        self.file_label.pack(side="left", padx=10)
        self.info_label = tk.Label(top, text="", **label_style)
        self.info_label.pack(side="right", padx=10)

        # Canvas for video
        self.canvas = tk.Canvas(self.root, bg="#000000", cursor="crosshair")
        self.canvas.pack(fill="both", expand=True, padx=5)
        self.canvas.bind("<ButtonPress-1>", self.on_mouse_down)
        self.canvas.bind("<B1-Motion>", self.on_mouse_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_mouse_up)
        self.canvas.bind("<ButtonPress-3>", self.clear_crop)

        # Timeline slider
        slider_frame = tk.Frame(self.root, **frame_style)
        slider_frame.pack(fill="x", padx=5)
        self.time_label = tk.Label(slider_frame, text="0:00.0 / 0:00.0", **label_style, font=("Consolas", 10))
        self.time_label.pack(side="left")
        self.slider = tk.Scale(slider_frame, from_=0, to=100, orient="horizontal",
                               command=self.on_slider, showvalue=False,
                               bg="#2b2b2b", fg="#ffffff", troughcolor="#555555",
                               highlightthickness=0, length=600)
        self.slider.pack(side="left", fill="x", expand=True, padx=5)

        # Playback controls
        ctrl = tk.Frame(self.root, **frame_style)
        ctrl.pack(fill="x", padx=5, pady=3)
        tk.Button(ctrl, text="<< 5s", command=lambda: self.seek_rel(-5), **btn_style).pack(side="left", padx=2)
        tk.Button(ctrl, text="< 1s", command=lambda: self.seek_rel(-1), **btn_style).pack(side="left", padx=2)
        tk.Button(ctrl, text="< Frame", command=lambda: self.seek_rel(-1/self.fps), **btn_style).pack(side="left", padx=2)
        self.play_btn = tk.Button(ctrl, text="Play", command=self.toggle_play, **btn_style, width=6)
        self.play_btn.pack(side="left", padx=6)
        tk.Button(ctrl, text="Frame >", command=lambda: self.seek_rel(1/self.fps), **btn_style).pack(side="left", padx=2)
        tk.Button(ctrl, text="1s >", command=lambda: self.seek_rel(1), **btn_style).pack(side="left", padx=2)
        tk.Button(ctrl, text="5s >>", command=lambda: self.seek_rel(5), **btn_style).pack(side="left", padx=2)

        # Trim/Crop/Speed controls
        edit_frame = tk.Frame(self.root, **frame_style)
        edit_frame.pack(fill="x", padx=5, pady=5)

        # Trim
        trim_lf = tk.LabelFrame(edit_frame, text="Trim", **label_style, font=("Segoe UI", 10, "bold"))
        trim_lf.pack(side="left", padx=5, fill="x", expand=True)
        row1 = tk.Frame(trim_lf, **frame_style)
        row1.pack(fill="x", padx=3, pady=2)
        tk.Label(row1, text="In:", **label_style).pack(side="left")
        self.in_var = tk.StringVar(value="0:00.000")
        self.in_entry = tk.Entry(row1, textvariable=self.in_var, width=10, bg="#333", fg="#fff",
                                 insertbackground="#fff")
        self.in_entry.pack(side="left", padx=3)
        tk.Button(row1, text="Set In", command=self.set_in, **btn_style).pack(side="left", padx=3)

        row2 = tk.Frame(trim_lf, **frame_style)
        row2.pack(fill="x", padx=3, pady=2)
        tk.Label(row2, text="Out:", **label_style).pack(side="left")
        self.out_var = tk.StringVar(value="0:00.000")
        self.out_entry = tk.Entry(row2, textvariable=self.out_var, width=10, bg="#333", fg="#fff",
                                  insertbackground="#fff")
        self.out_entry.pack(side="left", padx=3)
        tk.Button(row2, text="Set Out", command=self.set_out, **btn_style).pack(side="left", padx=3)
        tk.Button(row2, text="Reset Trim", command=self.reset_trim, **btn_style).pack(side="left", padx=3)

        # Crop
        crop_lf = tk.LabelFrame(edit_frame, text="Crop (drag on video, right-click=clear)", **label_style,
                                font=("Segoe UI", 10, "bold"))
        crop_lf.pack(side="left", padx=5, fill="x", expand=True)
        self.crop_label = tk.Label(crop_lf, text="No crop set", **label_style)
        self.crop_label.pack(padx=3, pady=2)

        # Speed
        speed_lf = tk.LabelFrame(edit_frame, text="Speed", **label_style, font=("Segoe UI", 10, "bold"))
        speed_lf.pack(side="left", padx=5)
        self.speed_var = tk.StringVar(value="1.0")
        tk.Entry(speed_lf, textvariable=self.speed_var, width=5, bg="#333", fg="#fff",
                 insertbackground="#fff").pack(padx=3, pady=2)
        tk.Label(speed_lf, text="x", **label_style).pack()

        # Export
        export_frame = tk.Frame(self.root, **frame_style)
        export_frame.pack(fill="x", padx=5, pady=5)
        tk.Button(export_frame, text="Export MP4", command=self.export, font=("Segoe UI", 11, "bold"),
                  bg="#2d7d2d", fg="#ffffff", activebackground="#3a9a3a", relief="flat",
                  padx=20, pady=6).pack(side="right", padx=5)
        self.status_label = tk.Label(export_frame, text="", **label_style)
        self.status_label.pack(side="left", padx=5)

    # --- File ---
    def open_file(self):
        path = filedialog.askopenfilename(filetypes=[
            ("Video", "*.mp4 *.avi *.mkv *.mov *.webm *.wmv"), ("All", "*.*")])
        if not path:
            return
        self.filepath = path
        if self.cap:
            self.cap.release()
        self.cap = cv2.VideoCapture(path)
        self.total_frames = int(self.cap.get(cv2.CAP_PROP_FRAME_COUNT))
        self.fps = self.cap.get(cv2.CAP_PROP_FPS) or 30
        self.width = int(self.cap.get(cv2.CAP_PROP_FRAME_WIDTH))
        self.height = int(self.cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
        self.duration = self.total_frames / self.fps

        self.file_label.config(text=os.path.basename(path))
        self.info_label.config(text=f"{self.width}x{self.height}  |  {self.fps:.1f} fps  |  {self.duration:.1f}s")

        self.slider.config(to=self.total_frames - 1)
        self.trim_in = 0.0
        self.trim_out = self.duration
        self.in_var.set(self._fmt_time(0))
        self.out_var.set(self._fmt_time(self.duration))
        self.speed_var.set("1.0")
        self.crop_rect = None
        self.crop_label.config(text="No crop set")
        self.current_frame = 0
        self.show_frame(0)

    # --- Display ---
    def show_frame(self, frame_num):
        if not self.cap:
            return
        frame_num = max(0, min(int(frame_num), self.total_frames - 1))
        self.current_frame = frame_num
        self.cap.set(cv2.CAP_PROP_POS_FRAMES, frame_num)
        ret, frame = self.cap.read()
        if not ret:
            return

        # Fit to canvas
        cw = self.canvas.winfo_width() or 800
        ch = self.canvas.winfo_height() or 500
        self.scale = min(cw / self.width, ch / self.height, 1.0)
        dw = int(self.width * self.scale)
        dh = int(self.height * self.scale)
        self.display_offset_x = (cw - dw) // 2
        self.display_offset_y = (ch - dh) // 2

        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        img = Image.fromarray(frame_rgb).resize((dw, dh), Image.LANCZOS)
        self.tk_img = ImageTk.PhotoImage(img)

        self.canvas.delete("all")
        self.canvas.create_image(self.display_offset_x, self.display_offset_y,
                                 anchor="nw", image=self.tk_img)

        # Draw crop rectangle
        if self.crop_rect:
            x1, y1, x2, y2 = self.crop_rect
            cx1 = int(x1 * self.scale) + self.display_offset_x
            cy1 = int(y1 * self.scale) + self.display_offset_y
            cx2 = int(x2 * self.scale) + self.display_offset_x
            cy2 = int(y2 * self.scale) + self.display_offset_y
            self.canvas.create_rectangle(cx1, cy1, cx2, cy2, outline="#00ff00", width=2, dash=(4, 4))

        # Draw trim markers on time label
        t = frame_num / self.fps
        self.time_label.config(text=f"{self._fmt_time(t)} / {self._fmt_time(self.duration)}")
        self.slider.set(frame_num)

    def _fmt_time(self, sec):
        m = int(sec) // 60
        s = sec - m * 60
        return f"{m}:{s:06.3f}"

    def _parse_time(self, text):
        """Parse M:SS.sss or SS.sss to seconds"""
        text = text.strip()
        if ":" in text:
            parts = text.split(":")
            return float(parts[0]) * 60 + float(parts[1])
        return float(text)

    # --- Slider / Seek ---
    def on_slider(self, val):
        if not self.playing:
            self.show_frame(int(float(val)))

    def seek_rel(self, dt):
        if not self.cap:
            return
        self.playing = False
        self.play_btn.config(text="Play")
        new_frame = self.current_frame + int(dt * self.fps)
        self.show_frame(new_frame)

    def toggle_play(self):
        if not self.cap:
            return
        self.playing = not self.playing
        self.play_btn.config(text="Pause" if self.playing else "Play")
        if self.playing:
            self._play_loop()

    def _play_loop(self):
        if not self.playing or not self.cap:
            return
        self.current_frame += 1
        if self.current_frame >= self.total_frames:
            self.playing = False
            self.play_btn.config(text="Play")
            return
        self.show_frame(self.current_frame)
        self.root.after(int(1000 / self.fps), self._play_loop)

    # --- Crop ---
    def _canvas_to_video(self, cx, cy):
        vx = (cx - self.display_offset_x) / self.scale
        vy = (cy - self.display_offset_y) / self.scale
        vx = max(0, min(vx, self.width))
        vy = max(0, min(vy, self.height))
        return vx, vy

    def on_mouse_down(self, event):
        if not self.cap:
            return
        self.dragging = True
        self.drag_start = (event.x, event.y)
        vx, vy = self._canvas_to_video(event.x, event.y)
        self.crop_start = (vx, vy)

    def on_mouse_drag(self, event):
        if not self.dragging:
            return
        vx1, vy1 = self.crop_start
        vx2, vy2 = self._canvas_to_video(event.x, event.y)
        self.crop_rect = (min(vx1, vx2), min(vy1, vy2), max(vx1, vx2), max(vy1, vy2))
        self.show_frame(self.current_frame)
        w = abs(vx2 - vx1)
        h = abs(vy2 - vy1)
        self.crop_label.config(text=f"{int(w)}x{int(h)} at ({int(min(vx1,vx2))},{int(min(vy1,vy2))})")

    def on_mouse_up(self, event):
        self.dragging = False
        if self.crop_rect:
            x1, y1, x2, y2 = self.crop_rect
            # Snap to even numbers (H.264 requirement)
            x1, y1 = int(x1) & ~1, int(y1) & ~1
            x2, y2 = (int(x2) + 1) & ~1, (int(y2) + 1) & ~1
            w, h = x2 - x1, y2 - y1
            if w < 10 or h < 10:
                self.crop_rect = None
                self.crop_label.config(text="No crop set")
            else:
                self.crop_rect = (x1, y1, x2, y2)
                self.crop_label.config(text=f"{w}x{h} at ({x1},{y1})")
            self.show_frame(self.current_frame)

    def clear_crop(self, event=None):
        self.crop_rect = None
        self.crop_label.config(text="No crop set")
        if self.cap:
            self.show_frame(self.current_frame)

    # --- Trim ---
    def set_in(self):
        if not self.cap:
            return
        t = self.current_frame / self.fps
        self.trim_in = t
        self.in_var.set(self._fmt_time(t))

    def set_out(self):
        if not self.cap:
            return
        t = self.current_frame / self.fps
        self.trim_out = t
        self.out_var.set(self._fmt_time(t))

    def reset_trim(self):
        self.trim_in = 0.0
        self.trim_out = self.duration
        self.in_var.set(self._fmt_time(0))
        self.out_var.set(self._fmt_time(self.duration))

    # --- Export ---
    def export(self):
        if not self.filepath:
            messagebox.showwarning("No file", "Open a video first.")
            return

        try:
            trim_in = self._parse_time(self.in_var.get())
            trim_out = self._parse_time(self.out_var.get())
            speed = float(self.speed_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid trim time or speed value.")
            return

        if speed <= 0 or speed > 100:
            messagebox.showerror("Error", "Speed must be between 0 and 100.")
            return

        if trim_out <= trim_in:
            messagebox.showerror("Error", f"Out ({trim_out:.3f}s) must be after In ({trim_in:.3f}s).")
            return

        print(f"Trim: {trim_in:.3f}s -> {trim_out:.3f}s  (duration: {trim_out - trim_in:.3f}s)  speed: {speed}x")

        base, ext = os.path.splitext(self.filepath)
        outpath = filedialog.asksaveasfilename(
            initialfile=os.path.basename(base) + "_edited.mp4",
            initialdir=os.path.dirname(self.filepath),
            defaultextension=".mp4",
            filetypes=[("MP4", "*.mp4")])
        if not outpath:
            return

        # Build FFmpeg command
        vfilters = []
        if self.crop_rect:
            x1, y1, x2, y2 = self.crop_rect
            w, h = int(x2 - x1), int(y2 - y1)
            vfilters.append(f"crop={w}:{h}:{int(x1)}:{int(y1)}")

        if speed != 1.0:
            vfilters.append(f"setpts={1.0/speed}*PTS")

        # Find ffmpeg binary
        ffmpeg_bin = shutil.which("ffmpeg")
        if not ffmpeg_bin:
            # Winget default install location
            ffmpeg_bin = os.path.expandvars(
                r"%LOCALAPPDATA%\Microsoft\WinGet\Packages\Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe\ffmpeg-8.1-full_build\bin\ffmpeg.exe")
        if not os.path.isfile(ffmpeg_bin):
            messagebox.showerror("Error", "Cannot find ffmpeg.exe. Make sure it's installed and on PATH.")
            return

        duration = trim_out - trim_in
        cmd = [ffmpeg_bin, "-y",
               "-ss", f"{trim_in:.3f}",
               "-t", f"{duration:.3f}",
               "-i", self.filepath,
               "-avoid_negative_ts", "make_zero"]

        if vfilters:
            cmd += ["-vf", ",".join(vfilters)]

        # Audio: handle speed change with atempo (supports 0.5-2.0 range, chain for wider)
        if speed != 1.0:
            atempo_filters = []
            s = speed
            while s > 2.0:
                atempo_filters.append("atempo=2.0")
                s /= 2.0
            while s < 0.5:
                atempo_filters.append("atempo=0.5")
                s /= 0.5
            atempo_filters.append(f"atempo={s:.4f}")
            cmd += ["-af", ",".join(atempo_filters)]

        cmd += ["-c:v", "libx264", "-preset", "medium", "-crf", "18",
                "-c:a", "aac", "-b:a", "192k",
                "-movflags", "+faststart",
                outpath]

        self.status_label.config(text="Exporting...")
        self.root.update()

        cmd_str = " ".join(cmd)
        print(f"FFmpeg command:\n{cmd_str}")

        def run_export():
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                if result.returncode == 0:
                    self.root.after(0, lambda: self.status_label.config(
                        text=f"Done: {os.path.basename(outpath)}"))
                    self.root.after(0, lambda: messagebox.showinfo("Done", f"Saved to:\n{outpath}"))
                else:
                    err = result.stderr[-1000:] if result.stderr else "Unknown error"
                    print(f"FFmpeg stderr:\n{err}")
                    self.root.after(0, lambda: self.status_label.config(text="Export failed"))
                    self.root.after(0, lambda e=err: messagebox.showerror("FFmpeg Error", e))
            except Exception as e:
                self.root.after(0, lambda: self.status_label.config(text="Export failed"))
                self.root.after(0, lambda e=str(e): messagebox.showerror("Error", e))

        threading.Thread(target=run_export, daemon=True).start()


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1000x750")
    root.minsize(800, 600)
    app = VideoTool(root)
    root.mainloop()
