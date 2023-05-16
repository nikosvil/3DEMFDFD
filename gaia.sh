#!/bin/bash

#================================================================================================
# GAIAHPC initialization script
# Usage: ./gaia.sh <config.dat 
#================================================================================================
# Authors: N. D. Vilanakis
# Details: Applied Mathematics and Computers Lab, Technical University of Crete
# * A 3D frequency-domain electromagnetic solver employing a high order compact finite-difference scheme
# * N. D. Vilanakis, N. Economou, E. N. Mathioudakis, A. Vafidis
# * Computers and Geosciences
#================================================================================================

# Check if stdin was redirected
if [ -t 0 ]

 then
# No redirection - warn and ask for manual mode
  echo "Warning: No configuration file was used as input."
  echo "Please supply a suitable configuration file by appending <config.dat to the script."
  echo "Alternatively, use manual mode to supply the parameters."
  while true; do
   read -p "Do you wish to enter manual mode? (y/n) " -n1 yn
   case $yn in
        [yY] ) echo " ";
               echo "Manual mode: Please input the following parameters:";
               break;;
        [nN] ) echo " ";
              echo "Exiting...";
               exit;;
        *  ) echo "";;
   esac
  done
 echo "Enter Lx:"
 read lx
 echo "Enter Ly:"
 read ly
 echo "Enter Lz:"
 read lz
 echo "Enter Tx:"
 read tx
 echo "Enter Ty:"
 read ty
 echo "Enter Tz:"
 read tz
 echo "Enter gl:"
 read gl
 echo "Enter sigma0:"
 read sigma0
 echo "Enter sigma:"
 read sigma
 echo "Enter frequency (Hz):"
 read f
 echo "All necessary parameters have been supplied."


else


# Automatic mode
 read lx ly lz
 read tx ty tz
 read gl
 read sigma0
 read sigma
 read f
fi

exec 0>&1 1>&0
# Read grid size
grid=`grep -o nx=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`\
"x"`grep -o ny=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`\
"x"`grep -o nz=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`
echo "Grid size now is: " $grid
echo " "
   while true; do
   read -p "Press c to change grid size, q to quit, or any other key to continue: " -n1 gs
   case $gs in
        [cC] ) echo " ";
               echo "Enter nx";
               read nx;
               echo "Enter ny";
               read ny;
               echo "Enter nz";
               read nz;
               sed -i 's/nx=[0-9]\+,/nx='$nx',/' ./halfspace/main.f;
               sed -i 's/ny=[0-9]\+,/ny='$ny',/' ./halfspace/main.f;
               sed -i 's/nz=[0-9]\+,/nz='$nz',/' ./halfspace/main.f;
#              sed -i 's/ns=[0-9]\+;/ns='$nx';/' ./repo/create_EbackFiles_Halfspace.m;
               sed -i 's/nx=[0-9]\+;/nx='$nx';/' ./repo/create_EbackFiles_Halfspace.m;
               sed -i 's/ny=[0-9]\+;/ny='$ny';/' ./repo/create_EbackFiles_Halfspace.m;
               sed -i 's/nz=[0-9]\+;/nz='$nz';/' ./repo/create_EbackFiles_Halfspace.m;
grid=`grep -o nx=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`\
"x"`grep -o ny=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`\
"x"`grep -o nz=[0-9]*, ./halfspace/main.f | cut -c 4- | rev | cut -c 2- | rev`;
echo "Grid size now is: " $grid;
               break;;
        [qQ] ) echo " ";
               echo "Exiting...";
               exit;;
        *  )  echo " ";
              break;;
   esac
  done

# Start working processes
echo `date`": Running Matlab..."
cd repo
/usr/local/MATLAB/R2014b/bin/matlab -nodesktop -nosplash -r "run('create_EbackFiles_Halfspace.m');exit;" | tail -n +17
cd ..
echo `date`": Matlab finished."
echo `date`": Running main exec with the following params:"
echo "                  " $lx $ly $lz $tx $ty $tz $gl $sigma0 $sigma $f
cd halfspace
make
cd ..
./halfspace/exe
echo `date`": Main exec finished."
echo `date`": Creating directory structure..."
dirstr=./halfspace/Results/"$grid"/"$lx"_"$ly"_"$lz"/"$tx"_"$ty"_"$tz"_"$gl"/"${sigma0}"_"${sigma}"_"$f"/`date +%Y-%m-%d-%H:%M`
echo $dirstr
mkdir -p $dirstr
mkdir -p $dirstr/f2m
cd repo
python3 convert_files.py
/usr/local/MATLAB/R2014b/bin/matlab -nodesktop -nosplash -r "run('convert_file_type.m');exit;" | tail -n +17
cd ../
echo `date`": Moving files..."
mv ./temp/EX.txt $dirstr
mv ./temp/EY.txt $dirstr
mv ./temp/EZ.txt $dirstr
mv ./temp/EtotalX.txt $dirstr
mv ./temp/EtotalY.txt $dirstr
mv ./temp/ETX.mat $dirstr/f2m 
mv ./temp/ETY.mat $dirstr/f2m
mv ./temp/Ex.mat $dirstr/f2m 
mv ./temp/Ey.mat $dirstr/f2m
mv ./temp/TotalEX.mat $dirstr/f2m 
mv ./temp/TotalEY.mat $dirstr/f2m
mv ./temp/info.txt $dirstr
echo `date`": Cleaning up..."
rm ./temp/*
echo `date`": Finished."
