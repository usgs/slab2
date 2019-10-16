#!/bin/bash
# Kirstie L Haynie
# Date started: 09/19/2019
# Bash shell script to make slab models
# Running this for all slabs can take up to 8 hours
# If there are any issues check the 3 letter slab code, the database date, the TodaysDate used, that the slab2 anaconda enviornment is activated, and the filepaths

if [ "$#" -ne 4 ]; then
    echo " "
    echo "Usage: bash MakeSlabModels.sh DataBaseDate TodaysDate NoCores SlabName/All"
    echo "Usage: All will run all slab models"
    echo "Usage: Or use the 3 letter slab model name to just run that model"
    echo "Example: bash MakeSlabModels.sh 04-18 09.19.19 6 All"
    echo "Example: bash MakeSlabModels.sh 04-18 09.19.19 6 alu"
    echo "Run from the slab2code directory"
    echo " "
    exit 1
fi

date
pwd

dbd=$1
date=$2
ncore=$3
slab=$4

if [ $slab == "All" ]; then
# Aleutian-Alaska:
echo "Working on the Aleutian-Alaska region..................."
python s2d.py -p alu -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f alu_${dbd}_input.csv
mv alu_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/aluinput.par -c ${ncore}
perl generalplot.pl alu ${date} Output/alu_slab2_${date} Input/alu_${dbd}_input.csv
perl generalxsec.pl alu ${date} Output/alu_slab2_${date} Input/alu_${dbd}_input.csv library/slab1grids/alu_slab1.0_clip.grd 50

# Calabria:
echo "Working on the Calabria region..................."
python s2d.py -p cal -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f cal_${dbd}_input.csv
mv cal_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/calinput.par -c ${ncore}
perl generalplot.pl cal ${date} Output/cal_slab2_${date} Input/cal_${dbd}_input.csv
perl generalxsec.pl cal ${date} Output/cal_slab2_${date} Input/cal_${dbd}_input.csv library/slab1grids/cal_slab1.0_clip.grd 50


# Caribbean:
echo "Working on the Caribbean region..................."
python s2d.py -p car -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f car_${dbd}_input.csv
mv car_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/carinput.par -c ${ncore}
perl generalplot.pl car ${date} Output/car_slab2_${date} Input/car_${dbd}_input.csv
perl generalxsec.pl car ${date} Output/car_slab2_${date} Input/car_${dbd}_input.csv library/slab1grids/car_slab1.0_clip.grd 50

# Central America:
echo "Working on the Central Amreica region..................."
python s2d.py -p cam -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f cam_${dbd}_input.csv
mv cam_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/caminput.par -c ${ncore}
perl generalplot.pl cam ${date} Output/cam_slab2_${date} Input/cam_${dbd}_input.csv/
perl generalxsec.pl cam ${date} Output/cam_slab2_${date} Input/cam_${dbd}_input.csv library/slab1grids/cam_slab1.0_clip.grd 50

# Cotabato:
echo "Working on the Cotabato region..................."
python s2d.py -p cot -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f cot_${dbd}_input.csv
mv cot_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/cotinput.par -c ${ncore}
perl generalplot.pl cot ${date} Output/cot_slab2_${date} Input/cot_${dbd}_input.csv/
perl generalxsec.pl cot ${date} Output/cot_slab2_${date} Input/cot_${dbd}_input.csv library/slab1grids/cot_slab1.0_clip.grd 50

# Cascadia:
echo "Working on the Cascadia region..................."
python s2d.py -p cas -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f cas_${dbd}_input.csv
mv cas_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/casinput.par -c ${ncore}
perl generalplot.pl cas ${date} Output/cas_slab2_${date} Input/cas_${dbd}_input.csv
perl generalxsec.pl cas ${date} Output/cas_slab2_${date} Input/cas_${dbd}_input.csv library/slab1grids/cas_slab1.0_clip.grd 50

# Halmahera:
echo "Working on the Halmahera region..................."
python s2d.py -p hal -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f hal_${dbd}_input.csv
mv hal_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/halinput.par -c ${ncore}
perl generalplot.pl hal ${date} Output/hal_slab2_${date} Input/hal_${dbd}_input.csv
perl generalxsec.pl hal ${date} Output/hal_slab2_${date} Input/hal_${dbd}_input.csv library/slab1grids/hal_slab1.0_clip.grd 50

# Hellenic:
echo "Working on the Hellenic region..................."
python s2d.py -p hel -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f hel_${dbd}_input.csv
mv hel_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/helinput.par -c ${ncore}
perl generalplot.pl hel ${date} Output/hel_slab2_${date} Input/hel_${dbd}_input.csv
perl generalxsec.pl hel ${date} Output/hel_slab2_${date} Input/hel_${dbd}_input.csv library/slab1grids/hel_slab1.0_clip.grd 50

# Himalaya:
echo "Working on the Himalaya ubduction zone..................."
python s2d.py -p him -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f him_${dbd}_input.csv
mv him_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/himinput.par -c ${ncore}
perl generalplot.pl him ${date} Output/him_slab2_${date} Input/him_${dbd}_input.csv
perl generalxsec.pl him ${date} Output/him_slab2_${date} Input/him_${dbd}_input.csv library/slab1grids/him_slab1.0_clip.grd 50

# Hindu Kush:
echo "Working on the Hindu Kush region..................."
python s2d.py -p hin -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f hin_${dbd}_input.csv
mv hin_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/hininput.par -c ${ncore}
perl generalplot.pl hin ${date} Output/hin_slab2_${date} Input/hin_${dbd}_input.csv
perl generalxsec.pl hin ${date} Output/hin_slab2_${date} Input/hin_${dbd}_input.csv library/slab1grids/hin_slab1.0_clip.grd 50

# Izu-Bonin:
echo "Working on the Izu-Bonin region..................."
python s2d.py -p izu -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f izu_${dbd}_input.csv
mv izu_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/izuinput.par -c ${ncore}
perl generalplot.pl izu ${date} Output/izu_slab2_${date} Input/izu_${dbd}_input.csv
perl generalxsec.pl izu ${date} Output/izu_slab2_${date} Input/izu_${dbd}_input.csv library/slab1grids/izu_slab1.0_clip.grd 50

# Kamchatka-Japng:
echo "Working on the Kamchatka-Japan region..................."
python s2d.py -p jap -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f jap_${dbd}_input.csv
mv jap_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/japinput.par -c ${ncore}
perl generalplot.pl jap ${date} Output/jap_slab2_${date} Input/jap_${dbd}_input.csv
perl generalxsec.pl jap ${date} Output/jap_slab2_${date} Input/jap_${dbd}_input.csv library/slab1grids/jap_slab1.0_clip.grd 50

# Kermadec:
echo "Working on the Kermadec region..................."
python s2d.py -p ker -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f ker_${dbd}_input.csv
mv ker_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/kerinput.par -c ${ncore}
perl generalplot.pl ker ${date} Output/ker_slab2_${date} Input/ker_${dbd}_input.csv
perl generalxsec.pl ker ${date} Output/ker_slab2_${date} Input/ker_${dbd}_input.csv library/slab1grids/ker_slab1.0_clip.grd 50

# Kuril:
echo "Working on the Kuril region..................."
python s2d.py -p kur -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f kur_${dbd}_input.csv
mv kur_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/kurinput.par -c ${ncore}
perl generalplot.pl kur ${date} Output/kur_slab2_${date} Input/kur_${dbd}_input.csv
perl generalxsec.pl kur ${date} Output/kur_slab2_${date} Input/kur_${dbd}_input.csv library/slab1grids/kur_slab1.0_clip.grd 50

# Makran:
echo "Working on the Makran region..................."
python s2d.py -p mak -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f mak_${dbd}_input.csv
mv mak_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/makinput.par -c ${ncore}
perl generalplot.pl mak ${date} Output/mak_slab2_${date} Input/mak_${dbd}_input.csv
perl generalxsec.pl mak ${date} Output/mak_slab2_${date} Input/mak_${dbd}_input.csv library/slab1grids/mak_slab1.0_clip.grd 50

# Manila:
echo "Working on the Manila Trench region..................."
python s2d.py -p man -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f man_${dbd}_input.csv
mv man_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/maninput.par -c ${ncore}
perl generalplot.pl man ${date} Output/man_slab2_${date} Input/man_${dbd}_input.csv
perl generalxsec.pl man ${date} Output/man_slab2_${date} Input/man_${dbd}_input.csv library/slab1grids/man_slab1.0_clip.grd 50

# Pamir:
echo "Working on the Pamir region..................."
python s2d.py -p pam -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f pam_${dbd}_input.csv
mv pam_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/paminput.par -c ${ncore}
perl generalplot.pl pam ${date} Output/pam_slab2_${date} Input/pam_${dbd}_input.csv
perl generalxsec.pl pam ${date} Output/pam_slab2_${date} Input/pam_${dbd}_input.csv library/slab1grids/pam_slab1.0_clip.grd 50

# Papa New Guinea:
echo "Working on the New Guinea region..................."
python s2d.py -p png -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f png_${dbd}_input.csv
mv png_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/pnginput.par -c ${ncore}
perl generalplot.pl png ${date} Output/png_slab2_${date} Input/png_${dbd}_input.csv
perl generalxsec.pl png ${date} Output/png_slab2_${date} Input/png_${dbd}_input.csv library/slab1grids/png_slab1.0_clip.grd 50

# Puysegur:
echo "Working on the Puysegur region..................."
python s2d.py -p puy -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f puy_${dbd}_input.csv
mv puy_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/puyinput.par -c ${ncore}
perl generalplot.pl puy ${date} Output/puy_slab2_${date} Input/puy_${dbd}_input.csv
perl generalxsec.pl puy ${date} Output/puy_slab2_${date} Input/puy_${dbd}_input.csv library/slab1grids/puy_slab1.0_clip.grd 50

# South American:
echo "Working on the South America region..................."
python s2d.py -p sam -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f sam_${dbd}_input.csv
mv sam_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/saminput.par -c ${ncore}
perl generalplot.pl sam ${date} Output/sam_slab2_${date} Input/sam_${dbd}_input.csv
perl generalxsec.pl sam ${date} Output/sam_slab2_${date} Input/sam_${dbd}_input.csv library/slab1grids/sam_slab1.0_clip.grd 50

# Scotia:
echo "Working on the Scotia Sea region..................."
python s2d.py -p sco -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f sco_${dbd}_input.csv
mv sco_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/scoinput.par -c ${ncore}
perl generalplot.pl sco ${date} Output/sco_slab2_${date} Input/sco_${dbd}_input.csv
perl generalxsec.pl sco ${date} Output/scou_slab2_${date} Input/sco_${dbd}_input.csv library/slab1grids/sco_slab1.0_clip.grd 50

# Solomon Islands:
echo "Working on the Solomon Islands region..................."
python s2d.py -p sol -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f sol_${dbd}_input.csv
mv sol_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/solinput.par -c ${ncore}
perl generalplot.pl sol ${date} Output/sol_slab2_${date} Input/sol_${dbd}_input.csv
perl generalxsec.pl sol ${date} Output/sol_slab2_${date} Input/sol_${dbd}_input.csv library/slab1grids/sol_slab1.0_clip.grd 50

# Sumatra-Java:
echo "Working on theSumatra-Java region..................."
python s2d.py -p sum -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f sum_${dbd}_input.csv
mv sum_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/suminput.par -c ${ncore}
perl generalplot.pl sum ${date} Output/sum_slab2_${date} Input/sum_${dbd}_input.csv/
perl generalxsec.pl sum ${date} Output/sum_slab2_${date} Input/sum_${dbd}_input.csv library/slab1grids/sum_slab1.0_clip.grd 50

# Vanuatu:
echo "Working on the Vanuatu region..................."
python s2d.py -p van -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f van_${dbd}_input.csv
mv van_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/vaninput.par -c ${ncore}
perl generalplot.pl van ${date} Output/van_slab2_${date} Input/van_${dbd}_input.csv
perl generalxsec.pl van ${date} Output/van_slab2_${date} Input/van_${dbd}_input.csv library/slab1grids/van_slab1.0_clip.grd 50

# Cotabato:
echo "Working on the Cotabato region..................."
python s2d.py -p cot -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f cot_${dbd}_input.csv
mv cot_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/cotinput.par -c ${ncore} -u Output/hal_slab2_${date}/hal_slab2_dep_${date}.grd
perl generalplot.pl cot ${date} Output/cot_slab2_${date} Input/cot_${dbd}_input.csv
perl generalxsec.pl cot ${date} Output/cot_slab2_${date} Input/cot_${dbd}_input.csv library/slab1grids/cot_slab1.0_clip.grd 50

# Muertos:
echo "Working on the Muertos Trough region..................."
python s2d.py -p mue -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f mue_${dbd}_input.csv
mv mue_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/mueinput.par -c ${ncore} -u Output/car_slab2_${date}/car_slab2_dep_${date}.grd
perl generalplot.pl mue ${date} Output/mue_slab2_${date} Input/mue_${dbd}_input.csv
perl generalxsec.pl mue ${date} Output/mue_slab2_${date} Input/mue_${dbd}_input.csv library/slab1gridsmue_slab1.0_clip.grd 50

# Sulawesi:
echo "Working on the Sulawesi region..................."
python s2d.py -p sul -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f sul_${dbd}_input.csv
mv sul_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/sulinput.par -c ${ncore} -u Output/hal_slab2_${date}/hal_slab2_dep_${date}.grd
perl generalplot.pl sul ${date} Output/sul_slab2_${date} Input/sul_${dbd}_input.csv
perl generalxsec.pl sul ${date} Output/sul_slab2_${date} Input/sul_${dbd}_input.csv library/slab1grids/sul_slab1.0_clip.grd 50

# Philippines:
echo "Working on the Philippines region..................."
python s2d.py -p phi -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f phi_${dbd}_input.csv
mv phi_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/phiinput.par -c ${ncore} -u Output/hal_slab2_${date}/hal_slab2_dep_${date}.grd
perl generalplot.pl phi ${date} Output/phi_slab2_${date} Input/phi_${dbd}_input.csv
perl generalxsec.pl phi ${date} Output/phi_slab2_${date} Input/phi_${dbd}_input.csv library/slab1grids/phi_slab1.0_clip.grd 50

# Ryukyu:
echo "Working on the Ryukyu region..................."
python s2d.py -p ryu -d /Users/khaynie/Src/Slab2/slab2-master/${dbd}database -f ryu_${dbd}_input.csv
mv ryu_${dbd}_input.csv Input
python slab2.py -p /Users/khaynie/Src/Slab2/slab2-master/slab2code/library/parameterfiles/ryuinput.par -c ${ncore} -u Output/kur_slab2_${date}/kur_slab2_dep_${date}.grd
perl generalplot.pl ryu ${date} Output/ryu_slab2_${date} Input/ryu_${dbd}_input.csv
perl generalxsec.pl ryu ${date} Output/ryu_slab2_${date} Input/ryu_${dbd}_input.csv library/slab1grids/ryu_slab1.0_clip.grd 50

echo "Finished with all slab models!!"
echo "To organize data into new directories and calculate the seismogenic zone thickness run:"
echo "python sztcalc.py -s Output -f Output/SZT-09.20.19 -i Input -t 04-18 -l o -d c -m 20 -x 65 -b e"
echo "FIRST: Change the file names in filelist on line 56 of sztcalc.py"
date
pwd

exit 1

elif [ $slab == "cot" ]; then
# specific slab
echo "Working on the ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore} -u Output/hal_slab2_09.19.19/hal_slab2_dep_09.19.19.grd
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with all slab models!!"
date
pwd

exit 1

elif [ $slab == "mue" ]; then
# specific slab
echo "Working on ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore} -u Output/car_slab2_09.19.19/car_slab2_dep_09.19.19.grd
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with the ${slab} model!!"
date
pwd

exit 1

elif [ $slab == "sul" ]; then
# specific slab
echo "Working on ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore} -u Output/hal_slab2_09.19.19/hal_slab2_dep_09.19.19.grd
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with the ${slab} model!!"
date
pwd

exit 1

elif [ $slab == "phi" ]; then
# specific slab
echo "Working on ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore} -u Output/hal_slab2_09.19.19/hal_slab2_dep_09.19.19.grd
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with the ${slab} model!!"
date
pwd

exit 1

elif [ $slab == "ryu" ]; then
# specific slab
echo "Working on ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore} -u Output/kur_slab2_09.19.19/kur_slab2_dep_09.19.19.grd
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with the ${slab} model!!"
date
pwd

exit 1

else
# specific slab
echo "Working on ${slab}..................."
python s2d.py -p ${slab} -d ../${dbd}database -f ${slab}_${dbd}_input.csv
mv ${slab}_${dbd}_input.csv Input
python slab2.py -p library/parameterfiles/${slab}input.par -c ${ncore}
perl generalplot.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv
perl generalxsec.pl ${slab} ${date} Output/${slab}_slab2_${date} Input/${slab}_${dbd}_input.csv library/slab1grids/${slab}_slab1.0_clip.grd 50

echo "Finished with the ${slab} model!!"
date
pwd

exit 1

fi
