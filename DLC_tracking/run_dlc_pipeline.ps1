# ===== USER SETTINGS =====
$config   = "C:\Users\Admin\Desktop\whisker_dorsal-RZ-2026-04-08\config.yaml"
$videoDir = "C:\260407_KA_electro_dorsal_whisking"
$python   = "C:\Users\Admin\.conda\envs\dlc310\python.exe"
$maxiters = 50000

# ===== RUN =====
& $python -c @"
import deeplabcut, glob, os

config = r'$config'
video_dir = r'$videoDir'

print('=== Creating training dataset ===')
deeplabcut.create_training_dataset(config)

print('=== Training ($maxiters iterations) ===')
deeplabcut.train_network(config, maxiters=$maxiters, saveiters=10000, displayiters=1000)

print('=== Evaluating ===')
deeplabcut.evaluate_network(config, plotting=True)

print('=== Analyzing videos ===')
vids = sorted(glob.glob(os.path.join(video_dir, '*.avi')))
print(f'Found {len(vids)} videos')
deeplabcut.analyze_videos(config, vids, save_as_csv=True)

print('\n=== Done ===')
print(f'CSV files saved in: {video_dir}')
"@

Read-Host "Press Enter to close"
