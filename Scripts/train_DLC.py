import deeplabcut;

camera = "Arm1";
monkey = "Skipper";
videoPath = "S:\Lab\Skipper\All Data\Skipper_02_28_2020\Videos\Camera_Arm 1\Renamed";

deeplabcut.create_new_project(camera,monkey,videoPath,copy_videos = False);
print("done");

