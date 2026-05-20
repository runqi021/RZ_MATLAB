"""DLC pipeline helper — called from MATLAB via system().

Usage:
    python dlc_pipeline.py create   --videoDir ... --projName ... --bodyparts vL1,vL2,...
    python dlc_pipeline.py extract  --config ...
    python dlc_pipeline.py label    --config ...
    python dlc_pipeline.py check    --config ...
    python dlc_pipeline.py train    --config ... [--net resnet_101] [--maxiters 200000]
    python dlc_pipeline.py evaluate --config ...
    python dlc_pipeline.py analyze  --config ... --videoDir ...
    python dlc_pipeline.py labeled  --config ... --videoDir ...
    python dlc_pipeline.py refine   --config ... --videoDir ...

    python dlc_pipeline.py full     --videoDir ... --projName ... --bodyparts vL1,vL2,...
        (runs create → extract → opens labeling GUI, then stops — resume with train after labeling)
"""

import argparse
import os
import glob
import deeplabcut


def find_videos(video_dir, exts=('avi', 'mp4', 'mov')):
    """Return list of video files in a directory."""
    vids = []
    for ext in exts:
        vids.extend(glob.glob(os.path.join(video_dir, f'*.{ext}')))
    vids.sort()
    return vids


def cmd_create(args):
    """Create a new DLC project and add videos."""
    vids = find_videos(args.videoDir)
    if not vids:
        print(f"ERROR: no videos found in {args.videoDir}")
        return None

    print(f"Found {len(vids)} videos in {args.videoDir}")

    bodyparts = [b.strip() for b in args.bodyparts.split(',')]
    print(f"Bodyparts ({len(bodyparts)}): {bodyparts}")

    # Pass only 1 video to create_new_project to avoid mass copying
    # (Windows lacks symlink privileges, so DLC copies every file)
    config_path = deeplabcut.create_new_project(
        args.projName,
        args.experimenter,
        [vids[0]],
        working_directory=args.workDir,
        copy_videos=False,
        multianimal=False,
    )
    print(f"\nProject created: {config_path}")

    # Manually add all videos to config (no copying) via config.yaml
    import ruamel.yaml as ry
    yaml_loader = ry.YAML()
    with open(config_path, 'r') as f:
        cfg_tmp = yaml_loader.load(f)
    # Point video_sets at the ORIGINAL directory — no copies needed
    cfg_tmp['video_sets'] = {}
    for v in vids:
        cfg_tmp['video_sets'][v] = {'crop': '0, 682, 0, 680'}
    with open(config_path, 'w') as f:
        yaml_loader.dump(cfg_tmp, f)
    # Remove the single copied/symlinked video in project/videos/
    proj_vid_dir = os.path.join(os.path.dirname(config_path), 'videos')
    for f in os.listdir(proj_vid_dir):
        fp = os.path.join(proj_vid_dir, f)
        if os.path.isfile(fp):
            os.remove(fp)
    print(f"Config updated: {len(vids)} videos added (no copies)")

    # Update config with bodyparts and skeleton
    import ruamel.yaml
    yaml = ruamel.yaml.YAML()
    with open(config_path, 'r') as f:
        cfg = yaml.load(f)

    cfg['bodyparts'] = bodyparts
    cfg['numframes2pick'] = args.numframes

    # Set network type
    cfg['net_type'] = args.net if hasattr(args, 'net') and args.net else 'resnet_101'

    with open(config_path, 'w') as f:
        yaml.dump(cfg, f)

    print(f"Config updated: bodyparts={bodyparts}, numframes2pick={args.numframes}")
    print(f"\nconfig_path = {config_path}")
    return config_path


def cmd_extract(args):
    """Extract frames using kmeans."""
    print(f"Extracting frames (algo={args.algo}, numframes2pick from config)...")
    deeplabcut.extract_frames(
        args.config,
        mode='automatic',
        algo=args.algo,
        userfeedback=False,
    )
    print("Frame extraction done.")


def cmd_label(args):
    """Open the labeling GUI."""
    print("Opening labeling GUI...")
    deeplabcut.label_frames(args.config)
    print("Labeling GUI closed.")


def cmd_check(args):
    """Check labels."""
    print("Checking labels...")
    deeplabcut.check_labels(args.config)
    print("Check done — see labeled-data folders for verification images.")


def cmd_train(args):
    """Create training dataset and train."""
    print("Creating training dataset...")
    deeplabcut.create_training_dataset(args.config)

    print(f"Training (maxiters={args.maxiters}, net={args.net})...")
    deeplabcut.train_network(
        args.config,
        maxiters=args.maxiters,
        saveiters=args.saveiters,
        displayiters=500,
    )
    print("Training done.")


def cmd_evaluate(args):
    """Evaluate the trained network."""
    print("Evaluating network...")
    deeplabcut.evaluate_network(args.config, plotting=True)
    print("Evaluation done.")


def cmd_analyze(args):
    """Run inference on videos."""
    vids = find_videos(args.videoDir)
    print(f"Analyzing {len(vids)} videos...")
    deeplabcut.analyze_videos(
        args.config,
        vids,
        save_as_csv=True,
    )
    print("Analysis done — CSV files saved alongside videos.")


def cmd_labeled(args):
    """Create labeled videos for QC."""
    vids = find_videos(args.videoDir)
    print(f"Creating labeled videos for {len(vids)} files...")
    deeplabcut.create_labeled_video(
        args.config,
        vids,
        draw_skeleton=True,
    )
    print("Labeled videos created.")


def cmd_refine(args):
    """Extract outlier frames for refinement."""
    vids = find_videos(args.videoDir)
    print("Extracting outlier frames...")
    deeplabcut.extract_outlier_frames(args.config, vids)
    print("Now run: python dlc_pipeline.py label --config ...")
    print("Then retrain.")


def cmd_full(args):
    """Create project, extract frames, open labeling GUI."""
    config_path = cmd_create(args)
    if config_path is None:
        return

    # Extract frames
    args.config = config_path
    args.algo = args.algo if hasattr(args, 'algo') else 'kmeans'
    cmd_extract(args)

    print("\n" + "="*60)
    print("NEXT STEPS:")
    print(f"  1. Label frames in the GUI that will open now")
    print(f"  2. Close GUI when done labeling")
    print(f"  3. Then run training from MATLAB:")
    print(f"     run_dlc('train', '{config_path}')")
    print("="*60 + "\n")

    cmd_label(args)


def main():
    parser = argparse.ArgumentParser(description="DLC pipeline helper")
    sub = parser.add_subparsers(dest='cmd')

    # --- create ---
    p = sub.add_parser('create')
    p.add_argument('--videoDir', required=True)
    p.add_argument('--projName', default='whisker_tracking')
    p.add_argument('--experimenter', default='RZ')
    p.add_argument('--bodyparts', required=True, help='comma-separated')
    p.add_argument('--workDir', default=None)
    p.add_argument('--numframes', type=int, default=20)
    p.add_argument('--net', default='resnet_101')

    # --- extract ---
    p = sub.add_parser('extract')
    p.add_argument('--config', required=True)
    p.add_argument('--algo', default='kmeans')

    # --- label ---
    p = sub.add_parser('label')
    p.add_argument('--config', required=True)

    # --- check ---
    p = sub.add_parser('check')
    p.add_argument('--config', required=True)

    # --- train ---
    p = sub.add_parser('train')
    p.add_argument('--config', required=True)
    p.add_argument('--net', default='resnet_101')
    p.add_argument('--maxiters', type=int, default=200000)
    p.add_argument('--saveiters', type=int, default=10000)

    # --- evaluate ---
    p = sub.add_parser('evaluate')
    p.add_argument('--config', required=True)

    # --- analyze ---
    p = sub.add_parser('analyze')
    p.add_argument('--config', required=True)
    p.add_argument('--videoDir', required=True)

    # --- labeled video ---
    p = sub.add_parser('labeled')
    p.add_argument('--config', required=True)
    p.add_argument('--videoDir', required=True)

    # --- refine ---
    p = sub.add_parser('refine')
    p.add_argument('--config', required=True)
    p.add_argument('--videoDir', required=True)

    # --- full (create+extract+label) ---
    p = sub.add_parser('full')
    p.add_argument('--videoDir', required=True)
    p.add_argument('--projName', default='whisker_tracking')
    p.add_argument('--experimenter', default='RZ')
    p.add_argument('--bodyparts', required=True)
    p.add_argument('--workDir', default=None)
    p.add_argument('--numframes', type=int, default=20)
    p.add_argument('--net', default='resnet_101')
    p.add_argument('--algo', default='kmeans')

    args = parser.parse_args()

    cmds = {
        'create': cmd_create,
        'extract': cmd_extract,
        'label': cmd_label,
        'check': cmd_check,
        'train': cmd_train,
        'evaluate': cmd_evaluate,
        'analyze': cmd_analyze,
        'labeled': cmd_labeled,
        'refine': cmd_refine,
        'full': cmd_full,
    }

    if args.cmd in cmds:
        cmds[args.cmd](args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
