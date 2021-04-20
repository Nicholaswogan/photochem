# delete stuff in build dir
cd ../cvode-5.7.0/build_dir
rm -r *
touch placeholder

# yaml dir
cd ../../fortran-yaml/build_dir
rm -r *
touch placeholder

cd ../../
rm -r StringiFor

cd dependencies
find . -type f ! -name '*.sh' -delete
rm -r -- ./*/
