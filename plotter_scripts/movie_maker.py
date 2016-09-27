
# A starting point?
#ffmpeg -y -i "frame%03d.png" movie_name.m4v
# framerate 1/5 means each image will last 5 seconds, 10 means 1/10 of a second.
# -i is for the input files
# -r is for the fps of the output movie
# -vframes is the number of images I wish to use.
# -start_number is the frame to start on
# -pattern_type glob allows usage of wildcard * if don't have all frames in sequence
# Not  101, 102 etc but rather 101, 103, 104 etc.
ffmpeg -framerate 10 -start_number 0170 -i movieframe_quad3_%04d_Density_2.0pc_65587.0.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4


exec("convert -geometry 1600x1600 -density 200x200 -quality 100 -resize 800x $pdf_path $temp_images");
exec("ffmpeg -loop 1 -f image2 -framerate 0.5 -b 1800 -i $temp_images_wildcard -c:v libx264 -preset slow -tune stillimage -r 5 -t 22 -y $frame_target 2>&1",$output);



ffmpeg -framerate 10 -start_number 179 -pattern_type glob -i '*65587*.png' -vframes 100 -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
