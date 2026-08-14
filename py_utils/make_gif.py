from pathlib import Path
from PIL import Image
import argparse

## This is called by animate.sh

def main():
    parser = argparse.ArgumentParser(description="Create a GIF from PNG frames.")
    parser.add_argument("frames_dir", type=str, help="Directory containing PNG frames")
    parser.add_argument("output_gif", type=str, help="Output GIF filename")
    parser.add_argument("--fps", type=float, default=2.0, help="Frames per second")
    args = parser.parse_args()

    frames_dir = Path(args.frames_dir)
    pngs = sorted(frames_dir.glob("*.png"))

    if not pngs:
        raise FileNotFoundError(f"No PNG frames found in {frames_dir}")

    duration_ms = int(1000 / args.fps)

    images = []
    base_size = None

    for png in pngs:
        im = Image.open(png).convert("RGBA")
        if base_size is None:
            base_size = im.size
        elif im.size != base_size:
            im = im.resize(base_size, Image.LANCZOS)
        images.append(im)

    first = images[0]
    rest = images[1:]

    first.save(
        args.output_gif,
        save_all=True,
        append_images=rest,
        duration=duration_ms,
        loop=0,
        disposal=2
    )

    print(f"Wrote GIF with {len(images)} frames to {args.output_gif}")

if __name__ == "__main__":
    main()
