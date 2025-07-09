### Pipeline Bash pour production et analyse/préparation de fichiers pour analyse, de simulation de dynamique moléculaire dans différentes conditions. 

#########################
#### Exemple d'index ####
#########################


gmx make_ndx -f maeva.gro -o index.ndx
###
1 | 13
name 18 SOLU

14
name 19 MEMB

17 | 15 | 16
name 20 SOLV

18 | 19

q

[ r_68_&_SG ]
1092
[ r_192_&_NE2 ]
3089
[ CU241 ]
3885
[ CU242 ]
3886
[ r_188_&_SG ]
3034
[ r_72_&_NE2 ]
1147

###
rm -f \#*\
####

################################
#### Preparation simulation ####
################################

#############
##### alone
#############
gmx grompp -f minimisation/step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -o minimisation/em_step5.tpr -maxwarn 1
gmx mdrun -v -deffnm minimisation/em_step5


gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_alone_1us.tpr
###prod test 
gmx grompp -f prod_test/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_test/md_test.tpr
gmx mdrun -v -deffnm prod_test/md_test.tpr



#############
##### ETHE
#############
gmx grompp -f minimisation/step6.0_minimization.mdp -c ETR2_ETHE.gro -r ETR2_ETHE.gro -p topol.top -o minimisation/em_step5.tpr
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###prod test 
gmx grompp -f prod_test/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_test/md_test.tpr
gmx mdrun -v -deffnm prod_test/md_test.tpr
###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_ETHE_1us.tpr


#############
##### ETOH
#############
gmx grompp -f minimisation/step6.0_minimization.mdp -c ETR2_ETOH.gro -r ETR2_ETOH.gro -p topol.top -o minimisation/em_step5.tpr
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###prod test 
gmx grompp -f prod_test/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_test/md_test.tpr
gmx mdrun -v -deffnm prod_test/md_test.tpr
###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_ETOH_1us.tpr


####################
#### ETHE & ETOH
####################
gmx grompp -f minimisation/step6.0_minimization.mdp -c ETR2_ETHE_ETOH.gro -r ETR2_ETHE_ETOH.gro -p topol.top -o minimisation/em_step5.tpr
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_BOTH_1us.tpr

###prod test 
gmx grompp -f prod_test/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_test/md_test.tpr
gmx mdrun -v -deffnm prod_test/md_test

###################################
######## Membrane Complexe ########
###################################

gmx editconf -f etoh.pdb -o ETOH.gro
gmx editconf -f ethe.pdb -o ETHE.gro

gmx insert-molecules -f step5_input.gro -ci ETHE.gro -replace TIP3 -nmol 10 -o ETR2_ETHE.gro
gmx insert-molecules -f ETR2_ETHE.gro -ci ETOH.gro -replace TIP3 -nmol 10 -o ETR2_ETHE_ETOH.gro


gmx grompp -f minimisation/step6.0_minimization.mdp -c ETR2_ETHE_ETOH.gro -r ETR2_ETHE_ETOH.gro -p topol.top -o minimisation/em_step5.tpr #-maxwarn 1
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/MEMB_complexe.tpr


###############################
######## 20 molécules  ########
###############################

gmx insert-molecules -f 20ETHE/step5_input.gro -ci 20ETHE/ETHE.gro -replace TIP3 -nmol 20 -o 20ETHE/20ETHE.gro
gmx insert-molecules -f 20ETOH/step5_input.gro -ci 20ETOH/ETOH.gro -replace TIP3 -nmol 20 -o 20ETOH/20ETOH.gro

gmx make_ndx -f 20ETOH/20ETOH.gro -o 20ETOH/index.ndx 
gmx make_ndx -f 20ETHE/20ETHE.gro -o 20ETHE/index.ndx 

#############
##### ETOH
#############
gmx grompp -f minimisation/step6.0_minimization.mdp -c 20ETOH.gro -r 20ETOH.gro -p topol.top -o minimisation/em_step5.tpr
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_20ETOH_1us.tpr


#############
##### ETHE
#############
gmx grompp -f minimisation/step6.0_minimization.mdp -c 20ETHE.gro -r 20ETHE.gro -p topol.top -o minimisation/em_step5.tpr
gmx mdrun -v -deffnm minimisation/em_step5

gmx grompp -f equil1/step6.1_equilibration.mdp -c minimisation/em_step5.gro -r minimisation/em_step5.gro -n index.ndx -p topol.top -o equil1/equil1.tpr
gmx mdrun -v -deffnm equil1/equil1 -nt 8 -pin on

gmx grompp -f equil2/step6.2_equilibration.mdp -c equil1/equil1.gro -r equil1/equil1.gro -n index.ndx -p topol.top -o equil2/equil2.tpr
gmx mdrun -v -deffnm equil2/equil2 -nt 8 -pin on

gmx grompp -f equil3/step6.3_equilibration.mdp -c equil2/equil2.gro -r equil2/equil2.gro -n index.ndx -p topol.top -o equil3/equil3.tpr
gmx mdrun -v -deffnm equil3/equil3 -nt 8 -pin on

gmx grompp -f equil4/step6.4_equilibration.mdp -c equil3/equil3.gro -r equil3/equil3.gro -n index.ndx -p topol.top -o equil4/equil4.tpr
gmx mdrun -v -deffnm equil4/equil4 -nt 8 -pin on

gmx grompp -f equil5/step6.5_equilibration.mdp -c equil4/equil4.gro -r equil4/equil4.gro -n index.ndx -p topol.top -o equil5/equil5.tpr
gmx mdrun -v -deffnm equil5/equil5 -nt 8 -pin on

gmx grompp -f equil6/step6.6_equilibration.mdp -c equil5/equil5.gro -r equil5/equil5.gro -n index.ndx -p topol.top -o equil6/equil6.tpr
gmx mdrun -v -deffnm equil6/equil6 -nt 8 -pin on

###1us
gmx grompp -f prod_1us/step7_production.mdp -c equil6/equil6.gro -r equil6/equil6.gro -n index.ndx -p topol.top -o prod_1us/md_20ETHE_1us.tpr




##########################################################################################################
#### Une fois que les simulations sont terminées et que nous pouvons travailler sur les productions : 
##########################################################################################################

###########################################
#### Supression conditions périodiques #### 
###########################################
#L'idée principale est de corriger les trajectoires pour que les molécules ou les complexes 
#ne soient pas artificiellement séparés à cause des bords de la boîte périodique.

# 1. Corriger les sauts d'abord
echo -e "0" | gmx trjconv -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us.xtc -o alone/md_VvETR2_1us_nojump.xtc -pbc nojump
echo -e "0" | gmx trjconv -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us.xtc -o ethylene/md_ETHE_1us_nojump.xtc -pbc nojump
echo -e "0" | gmx trjconv -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us.xtc -o ethanol/md_ETOH_1us_nojump.xtc -pbc nojump
echo -e "0" | gmx trjconv -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us.xtc -o both/md_ETHE_ETOH_1us_nojump.xtc -pbc nojump
# systeme

# 2. Centrer la protéine (très important AVANT le reemboîtage)
echo -e " 1" "0" | gmx trjconv -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_nojump.xtc -o alone/md_VvETR2_1us_center.xtc -center
echo -e " 1" "0" | gmx trjconv -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_nojump.xtc -o ethylene/md_ETHE_1us_center.xtc -center
echo -e " 1" "0" | gmx trjconv -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_nojump.xtc -o ethanol/md_ETOH_1us_center.xtc -center
echo -e " 1" "0" | gmx trjconv -s  both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_nojump.xtc -o both/md_ETHE_ETOH_1us_center.xtc -center


# 3. Réemboîter membrane + protéine
echo -e "0" | gmx trjconv -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_center.xtc -o alone/md_VvETR2_1us_final.xtc -pbc mol -ur compact
echo -e "0" | gmx trjconv -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_center.xtc -o ethylene/md_ETHE_1us_final.xtc -pbc mol -ur compact
echo -e "0" | gmx trjconv -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_center.xtc -o ethanol/md_ETOH_1us_final.xtc -pbc mol -ur compact
echo -e "0" | gmx trjconv -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_center.xtc -o both/md_ETHE_ETOH_1us_final.xtc -pbc mol -ur compact
#systeme

### TOUT TESTER ###
vmd alone/md_VvETR2_1us.gro alone/md_VvETR2_1us_final.xtc
vmd ethylene/md_ETHE_1us.gro ethylene/md_ETHE_1us_final.xtc
vmd ethanol/md_ETOH_1us.gro ethanol/md_ETOH_1us_final.xtc
vmd both/md_ETHE_ETOH_1us.gro both/md_ETHE_ETOH_1us_final.xtc




###############
####  Calcul RMSD : 
###############
gmx make_ndx -f  alone/md_VvETR2_1us.tpr -o alone/helices.ndx
gmx make_ndx -f  ethylene/md_ETHE_1us.tpr -o ethylene/helices.ndx
gmx make_ndx -f  ethanol/md_ETOH_1us.tpr -o ethanol/helices.ndx
gmx make_ndx -f  both/md_ETHE_ETOH_1us.tpr -o both/helices.ndx

r 18-49 | r 55-77 | r 84-118 | r 138-169 | r 175-197 | r 204-238
name 19 HELICES
q

#protein
echo -e "1" "1" | gmx rms -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -o alone/rmsd_VvETR2_prot.xvg
echo -e "1" "1" | gmx rms -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -o ethylene/rmsd_ETHE_prot.xvg
echo -e "1" "1" | gmx rms -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -o ethanol/rmsd_ETOH_prot.xvg
echo -e "1" "1" | gmx rms -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -o both/rmsd_both_prot.xvg


echo -e "3" "3" | gmx rms -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -n alone/helices.ndx -o alone/rmsd_ETR2Calpha.xvg
echo -e "3" "3" | gmx rms -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -n ethylene/helices.ndx -o ethylene/rmsd_ETHECalpha.xvg
echo -e "3" "3" | gmx rms -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/helices.ndx -o ethanol/rmsd_ETOHCalpha.xvg
echo -e "3" "3" | gmx rms -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -n both/helices.ndx -o both/rmsd_bothCalpha.xvg


echo -e "18" "18"| gmx rms -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -n alone/helices.ndx -o alone/rmsd_ETR2_helices.xvg
echo -e "19" "19"| gmx rms -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -n ethylene/helices.ndx -o ethylene/rmsd_ETHE_helices.xvg
echo -e "19" "19"| gmx rms -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/helices.ndx -o ethanol/rmsd_ETOH_helices.xvg
echo -e "20" "20"| gmx rms -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -n both/helices.ndx -o both/rmsd_both_helices.xvg

#mettre en nanosecondes : 

for file in alone/rmsd_VvETR2_prot.xvg ethylene/rmsd_ETHE_prot.xvg ethanol/rmsd_ETOH_prot.xvg both/rmsd_both_prot.xvg alone/rmsd_ETR2Calpha.xvg ethylene/rmsd_ETHECalpha.xvg ethanol/rmsd_ETOHCalpha.xvg both/rmsd_bothCalpha.xvg alone/rmsd_ETR2_helices.xvg ethylene/rmsd_ETHE_helices.xvg ethanol/rmsd_ETOH_helices.xvg both/rmsd_both_helices.xvg; do
    out="${file%.xvg}_ns.xvg"
    awk '
        $1 ~ /^@/ || $1 ~ /^#/ {
            gsub("ps", "ns")  # Remplace "ps" par "ns" 
            print
        }
        $1 !~ /^[@#]/ {
            printf "%.6f %s\n", $1/1000, $2
        }
    ' "$file" > "$out"
    echo "✔ Fichier converti : $out"
done


xmgrace alone/rmsd_VvETR2_prot_ns.xvg alone/rmsd_ETR2Calpha_ns.xvg alone/rmsd_ETR2_helices_ns.xvg
#alone
xmgrace ethylene/rmsd_ETHE_prot_ns.xvg ethylene/rmsd_ETHECalpha_ns.xvg ethylene/rmsd_ETHE_helices_ns.xvg
#ETHE
xmgrace ethanol/rmsd_ETOH_prot_ns.xvg ethanol/rmsd_ETOHCalpha_ns.xvg ethanol/rmsd_ETOH_helices_ns.xvg
#ETOH
xmgrace both/rmsd_both_prot_ns.xvg both/rmsd_bothCalpha_ns.xvg both/rmsd_both_helices_ns.xvg
#BOTH

xmgrace alone/rmsd_VvETR2_prot_ns.xvg ethylene/rmsd_ETHE_prot_ns.xvg ethanol/rmsd_ETOH_prot_ns.xvg both/rmsd_both_prot_ns.xvg
#RMSD evolution of all protein atoms of protein VvETR2 under different experimental conditions

xmgrace alone/rmsd_ETR2Calpha_ns.xvg ethylene/rmsd_ETHECalpha_ns.xvg ethanol/rmsd_ETOHCalpha_ns.xvg both/rmsd_bothCalpha_ns.xvg
#RMSD evolution of C-alpha atoms of protein VvETR2 across different experimental conditions 

xmgrace alone/rmsd_ETR2_helices_ns.xvg ethylene/rmsd_ETHE_helices_ns.xvg ethanol/rmsd_ETOH_helices_ns.xvg both/rmsd_both_helices_ns.xvg
#RMSD evolution of atoms in the helical region of protein VvETR2 under different experimental conditions


#####################
#### RMSD SMOOTH #### 
#####################
#Possible après traitement python
echo -e "3" "3" | gmx rms -s alone/md_VvETR2_1us.tpr -f alone/alone_smooth.xtc -n alone/helices.ndx -o alone/rmsd_ETR2Calpha_smooth.xvg
echo -e "3" "3" | gmx rms -s ethylene/md_ETHE_1us.tpr -f ethylene/ETHE_smooth.xtc -n ethylene/helices.ndx -o ethylene/rmsd_ETHECalpha_smooth.xvg
echo -e "3" "3" | gmx rms -s ethanol/md_ETOH_1us.tpr -f ethanol/ETOH_smooth.xtc -n ethanol/helices.ndx -o ethanol/rmsd_ETOHCalpha_smooth.xvg
echo -e "3" "3" | gmx rms -s both/md_ETHE_ETOH_1us.tpr -f both/BOTH_smooth.xtc -n both/helices.ndx -o both/rmsd_bothCalpha_smooth.xvg

echo -e "18" "18"| gmx rms -s alone/md_VvETR2_1us.tpr -f alone/alone_smooth.xtc -n alone/helices.ndx -o alone/rmsd_ETR2_helices_smooth.xvg
echo -e "19" "19"| gmx rms -s ethylene/md_ETHE_1us.tpr -f ethylene/ETHE_smooth.xtc -n ethylene/helices.ndx -o ethylene/rmsd_ETHE_helices_smooth.xvg
echo -e "19" "19"| gmx rms -s ethanol/md_ETOH_1us.tpr -f ethanol/ETOH_smooth.xtc -n ethanol/helices.ndx -o ethanol/rmsd_ETOH_helices_smooth.xvg
echo -e "20" "20"| gmx rms -s both/md_ETHE_ETOH_1us.tpr -f both/BOTH_smooth.xtc -n both/helices.ndx -o both/rmsd_both_helices_smooth.xvg

### METTRE EN ns
for file in alone/rmsd_ETR2Calpha_smooth.xvg ethylene/rmsd_ETHECalpha_smooth.xvg ethanol/rmsd_ETOHCalpha_smooth.xvg both/rmsd_bothCalpha_smooth.xvg alone/rmsd_ETR2_helices_smooth.xvg ethylene/rmsd_ETHE_helices_smooth.xvg ethanol/rmsd_ETOH_helices_smooth.xvg both/rmsd_both_helices_smooth.xvg; do
    out="${file%.xvg}_ns.xvg"
    awk '$1 ~ /^@/ || $1 ~ /^#/ {print} $1 !~ /^[@#]/ {printf "%.6f %s\n", $1/1000, $2}' "$file" > "$out"
    echo "✔ Fichier converti : $out"
done

xmgrace alone/rmsd_ETR2Calpha_smooth_ns.xvg ethylene/rmsd_ETHECalpha_smooth_ns.xvg ethanol/rmsd_ETOHCalpha_smooth_ns.xvg both/rmsd_bothCalpha_smooth_ns.xvg
#RMSD evolution of C-alpha atoms of protein VvETR2 across different experimental conditions, from smoothed trajectories

xmgrace alone/rmsd_ETR2_helices_smooth_ns.xvg ethylene/rmsd_ETHE_helices_smooth_ns.xvg ethanol/rmsd_ETOH_helices_smooth_ns.xvg both/rmsd_both_helices_smooth_ns.xvg
#RMSD evolution of atoms in the helical region of protein VvETR2 under different experimental conditions, from smoothed trajectories

###############
#### Calcul RMSF
###############
#carbone alpha
echo -e "3" "3" | gmx rmsf -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -o alone/rmsf_VvETR2_Calpha.xvg -res
echo -e "3" "3" | gmx rmsf -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -o ethylene/rmsf_ETHE_Calpha.xvg -res
echo -e "3" "3" | gmx rmsf -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -o ethanol/rmsf_ETOH_Calpha.xvg -res
echo -e "3" "3" | gmx rmsf -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -o both/rmsf_both_Calpha.xvg -res

#prot
echo -e "1" "1" | gmx rmsf -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -o alone/rmsf_VvETR2_prot.xvg -res
echo -e "1" "1" | gmx rmsf -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -o ethylene/rmsf_ETHE_prot.xvg -res
echo -e "1" "1" | gmx rmsf -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -o ethanol/rmsf_ETOH_prot.xvg -res
echo -e "1" "1" | gmx rmsf -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -o both/rmsf_both_prot.xvg -res

xmgrace alone/rmsf_VvETR2_prot.xvg ethylene/rmsf_ETHE_prot.xvg ethanol/rmsf_ETOH_prot.xvg both/rmsf_both_prot.xvg
#RMSF of all protein atoms of VvETR2 across different experimental conditions

#helices NON
echo -e "18" "18" | gmx rmsf -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_final.xtc -n alone/helices.ndx -o alone/rmsf_VvETR2_helix.xvg -res
echo -e "19" "19" | gmx rmsf -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -n ethylene/helices.ndx -o ethylene/rmsf_ETHE_helix.xvg -res
echo -e "19" "19" | gmx rmsf -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/helices.ndx -o ethanol/rmsf_ETOH_helix.xvg -res
echo -e "20" "20" | gmx rmsf -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -n both/helices.ndx -o both/rmsf_both_helix.xvg -res

xmgrace alone/rmsf_VvETR2_helix.xvg ethylene/rmsf_ETHE_helix.xvg ethanol/rmsf_ETOH_helix.xvg both/rmsf_both_helix.xvg

###############
###############
#### HOLE #### NE PAS FAIRE
gmx trjconv -f alone/md_VvETR2_1us_final.xtc -s alone/md_VvETR2_1us.tpr -o alone/hole/pore_frame0.pdb -dump 0
gmx trjconv -f alone/md_VvETR2_1us_final.xtc -s alone/md_VvETR2_1us.tpr -o alone/hole/pore_frame1.pdb -dump 1000
...

for i in {0..9}; do
    time_ns=$((i * 1000))  # temps en ps
    gmx trjconv -f alone/md_VvETR2_1us_final.xtc \
                -s alone/md_VvETR2_1us.tpr \
                -o alone/hole/pore_frame${i}.pdb \
                -dump ${time_ns} <<< "0"
done

#dans le bon dossier
for f in *.pdb; do
    grep '^ATOM' "$f" | grep -v 'TIP3\|DOPC\|CLA\|POT' > only_protein_"$f"
done

###SphereCon

spherecon.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--bb] [--dm] -i: Path to an input file in PDB (Protein Data Bank) file format or distance matrix file in case of --dm.

for f in only_protein*.pdb






##################
##################
#fit pour le hole
echo -e "3" "0" | gmx trjconv -fit rot+trans -f alone/md_VvETR2_1us_final.xtc -s alone/md_VvETR2_1us.tpr -o alone/md_VvETR2_1us_fit.xtc
echo -e "3" "0" | gmx trjconv -fit rot+trans -f ethylene/md_ETHE_1us_final.xtc -s ethylene/md_ETHE_1us.tpr -o ethylene/md_ETHE_1us_fit.xtc
echo -e "3" "0" | gmx trjconv -fit rot+trans -f ethanol/md_ETOH_1us_final.xtc -s ethanol/md_ETOH_1us.tpr -o ethanol/md_ETOH_1us_fit.xtc
echo -e "3" "0" | gmx trjconv -fit rot+trans -f both/md_ETHE_ETOH_1us_final.xtc -s both/md_ETHE_ETOH_1us.tpr -o both/md_ETHE_ETOH_1us_fit.xtc

#####################
#### DIFFUSION MB 
#####################

### Def boite proteine
# Sélection de la protéine sur vmd tout ça : 
set sel [atomselect top "protein"]
set minmax [measure minmax $sel]
set min [lindex $minmax 0]
set max [lindex $minmax 1]

# Extraire les coordonnées
set xmin [lindex $min 0]
set ymin [lindex $min 1]
set zmin [lindex $min 2]
set xmax [lindex $max 0]
set ymax [lindex $max 1]
set zmax [lindex $max 2]

# Définir les coins de la boîte
set p1 [list $xmin $ymin $zmin]
set p2 [list $xmax $ymin $zmin]
set p3 [list $xmin $ymax $zmin]
set p4 [list $xmax $ymax $zmin]
set p5 [list $xmin $ymin $zmax]
set p6 [list $xmax $ymin $zmax]
set p7 [list $xmin $ymax $zmax]
set p8 [list $xmax $ymax $zmax]

# Couleur rouge
draw color red

# Dessiner les 12 arêtes
draw line $p1 $p2
draw line $p1 $p3
draw line $p1 $p5

draw line $p2 $p4
draw line $p2 $p6

draw line $p3 $p4
draw line $p3 $p7

draw line $p4 $p8

draw line $p5 $p6
draw line $p5 $p7

draw line $p6 $p8
draw line $p7 $p8

##On réduit
# Sélection de la protéine
set sel [atomselect top "protein"]
set minmax [measure minmax $sel]
set min [lindex $minmax 0]
set max [lindex $minmax 1]

# Extraire les coordonnées
set xmin [lindex $min 0]
set ymin [lindex $min 1]
set zmin [lindex $min 2]
set xmax [lindex $max 0]
set ymax [lindex $max 1]
set zmax [lindex $max 2]

# Réduction : plus grand padding = boîte plus serrée
set padding 12.0

set xmin [expr $xmin + $padding]
set ymin [expr $ymin + $padding]
set zmin [expr $zmin + $padding]
set xmax [expr $xmax - $padding]
set ymax [expr $ymax - $padding]
set zmax [expr $zmax - $padding]

# Coins de la boîte
set p1 [list $xmin $ymin $zmin]
set p2 [list $xmax $ymin $zmin]
set p3 [list $xmin $ymax $zmin]
set p4 [list $xmax $ymax $zmin]
set p5 [list $xmin $ymin $zmax]
set p6 [list $xmax $ymin $zmax]
set p7 [list $xmin $ymax $zmax]
set p8 [list $xmax $ymax $zmax]

draw color pink

# Arêtes de la boîte
draw line $p1 $p2
draw line $p1 $p3
draw line $p1 $p5


draw line $p2 $p4
draw line $p2 $p6

draw line $p3 $p4
draw line $p3 $p7

draw line $p4 $p8

draw line $p5 $p6
draw line $p5 $p7

draw line $p6 $p8
draw line $p7 $p8

####
X : de 18.58 à 60.36

Y : de 21.40 à 57.54

Z : de 20.68 à 92.10





####

grep "ETOH" ethanol/md_ETOH_1us.gro 

#test ethanol 
gmx make_ndx -f ethanol/md_ETOH_1us.tpr -o ethanol/ethanol_only.ndx << EOF
r 15834
r 15835
r 15836
r 15837
r 15838
r 15839
r 15840
r 15841
r 15842
r 15843
q
EOF

for i in {19..28}; do
  outfile="ethanol/traj_etoh/ethanol_z$((i-18)).xvg"
  echo -e "$i" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox "$outfile" -com -ng 1
done

echo -e "19" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z1.xvg -com -ng 1
echo -e "20" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z2.xvg -com -ng 1
echo -e "21" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z3.xvg -com -ng 1
echo -e "22" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z4.xvg -com -ng 1
echo -e "23" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z5.xvg -com -ng 1
echo -e "24" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z6.xvg -com -ng 1
echo -e "25" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z7.xvg -com -ng 1
echo -e "26" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z8.xvg -com -ng 1
echo -e "27" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z9.xvg -com -ng 1
echo -e "28" | gmx traj -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_final.xtc -n ethanol/ethanol_only.ndx -ox ethanol/traj_etoh/ethanol_z10.xvg -com -ng 1

#apres les trajectoires, on applique le code python Suivre DeltaZ

#pour ethylène : 
grep "ETHE" ethylene/md_ETHE_1us.gro 

gmx make_ndx -f ethylene/md_ETHE_1us.tpr -o ethylene/ethylene_only.ndx << EOF
r 15850
r 15851
r 15852
r 15853
r 15854
r 15855
r 15856
r 15857
r 15858
r 15859
q
EOF

for i in {19..28}; do
  outfile="ethylene/traj_ethe/ethylene_z$((i-18)).xvg"
  echo -e "$i" | gmx traj -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_final.xtc -n ethylene/ethylene_only.ndx -ox "$outfile" -com -ng 1
done

###Pour les deux : 
grep "ETOH" both/md_ETHE_ETOH_1us.gro 
grep "ETHE" both/md_ETHE_ETOH_1us.gro 

###ajuster
gmx make_ndx -f both/md_ETHE_ETOH_1us.tpr -o both/index_ligands.ndx << EOF
r 15814
r 15815
r 15816
r 15817
r 15818
r 15819
r 15820
r 15821
r 15822
r 15823
r 15824
r 15825
r 15826
r 15827
r 15828
r 15829
r 15830
r 15831
r 15832
r 15833
q
EOF

###ajuster tout
for i in {20..29}; do
  outfile="both/traj_ethe/ethylene_z$((i-19)).xvg"
  echo -e "$i" | gmx traj -s  both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -n both/index_ligands.ndx -ox "$outfile" -com -ng 1
done

for i in {30..39}; do
  outfile="both/traj_etoh/ethanol_z$((i-29)).xvg"
  echo -e "$i" | gmx traj -s  both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_final.xtc -n both/index_ligands.ndx -ox "$outfile" -com -ng 1
done


###################
#### Distances ####
###################

#### ETHYLENE ####

# Fichiers d'entrée
XTC="ethylene/md_ETHE_1us_final.xtc"
TPR="ethylene/md_ETHE_1us.tpr"
SCRIPT="/home/mabadielim/Desktop/script/dist_Maeva.py"
OUTDIR="ethylene/Smooth_all_distances"
mkdir -p "$OUTDIR"

# Boucle sur les résidus (de 1 à 240)
for resid in $(seq 1 240); do
    OUTFILE="${OUTDIR}/resid_${resid}_allethylene.svg"
    echo "Processing resid $resid..."
    python3 "$SCRIPT" -f "$XTC" -s "$TPR" --select1 "resname ETHE" --select2 "resid $resid" -o "$OUTFILE"
done

#### ETHANOL ####

# Fichiers d'entrée
XTC="ethanol/md_ETOH_1us_final.xtc"
TPR="ethanol/md_ETOH_1us.tpr"
SCRIPT="/home/mabadielim/Desktop/script/dist_Maeva.py"
OUTDIR="ethanol/all_distances"
mkdir -p "$OUTDIR"

# Boucle sur les résidus (de 1 à 2)
for resid in $(seq 1 240); do
    OUTFILE="${OUTDIR}/resid_${resid}_allethanol.svg"
    echo "Processing resid $resid..."
    python3 "$SCRIPT" -f "$XTC" -s "$TPR" --select1 "resname ETOH" --select2 "resid $resid" -o "$OUTFILE"
done

#### BOTH ####

XTC="both/md_ETHE_ETOH_1us_final.xtc"
TPR="both/md_ETHE_ETOH_1us.tpr"
SCRIPT="/home/mabadielim/Desktop/script/dist_Maeva.py"
OUTDIR="both/Smooth_ETHE_all_distances"
mkdir -p "$OUTDIR"

# Boucle sur les résidus (de 1 à 240)
for resid in $(seq 1 240); do
    OUTFILE="${OUTDIR}/resid_${resid}_allethylene.svg"
    echo "Processing resid $resid..."
    python3 "$SCRIPT" -f "$XTC" -s "$TPR" --select1 "resname ETHE" --select2 "resid $resid" -o "$OUTFILE"
done

XTC="both/md_ETHE_ETOH_1us_final.xtc"
TPR="both/md_ETHE_ETOH_1us.tpr"
SCRIPT="/home/mabadielim/Desktop/script/dist_Maeva.py"
OUTDIR="both/Smooth_ETOH_all_distances"
mkdir -p "$OUTDIR"

# Boucle sur les résidus (de 1 à 240)
for resid in $(seq 1 240); do
    OUTFILE="${OUTDIR}/resid_${resid}_allethanol.svg"
    echo "Processing resid $resid..."
    python3 "$SCRIPT" -f "$XTC" -s "$TPR" --select1 "resname ETOH" --select2 "resid $resid" -o "$OUTFILE"
done


############################
#### Rayons de gyration ####
############################

#######
# ALONE
gmx make_ndx -f alone/md_VvETR2_1us.tpr -o alone/index_rg.ndx << EOF
r 23 - 28
r 143 - 148
r 88 - 93
r 208 - 213
r 74 - 79
r 194 - 199
20 | 21 | 22 | 23 | 18 | 19
name 24 HAUT
del 18
del 18 
del 18 
del 18 
del 18 
del 18 
r 33 - 38
r 153 - 158
r 65 - 70
r 185 - 190
r 99 - 104
r 219 - 224
20 | 21 | 22 | 23 | 24 | 19
del 19 
del 19 
del 19 
del 19 
del 19 
del 19 
name 19 MILIEU
r 43 - 48
r 163 - 168
r 55 - 60
r 175 - 180
r 110 - 115
r 230 - 235
20 | 21 | 22 | 23 | 24 | 25
del 20
del 20
del 20
del 20
del 20
del 20
name 20 BAS
r 18 - 49
name 21 H1
r 55 - 77 
name 22 H2
r 84-117
name 23 H3
r 138- 169
name 24 H4
r 175- 197
name 25 H5
r 204-237
name 26 H6
q
EOF

#######
# ETHE
gmx make_ndx -f ethylene/md_ETHE_1us.tpr -o ethylene/index_rg.ndx << EOF
r 23 - 28
r 143 - 148
r 88 - 93
r 208 - 213
r 74 - 79
r 194 - 199
20 | 21 | 22 | 23 | 24 | 19
name 25 HAUT
del 19 
del 19 
del 19 
del 19 
del 19 
del 19 
r 33 - 38
r 153 - 158
r 65 - 70
r 185 - 190
r 99 - 104
r 219 - 224
20 | 21 | 22 | 23 | 24 | 25
del 20
del 20
del 20
del 20
del 20
del 20
name 20 MILIEU
r 43 - 48
r 163 - 168
r 55 - 60
r 175 - 180
r 110 - 115
r 230 - 235
26 | 21 | 22 | 23 | 24 | 25
del 21
del 21
del 21
del 21
del 21
del 21
name 21 BAS
r 18 - 49
name 22 H1
r 55 - 77 
name 23 H2
r 84-117
name 24 H3
r 138- 169
name 25 H4
r 175- 197
name 26 H5
r 204-237
name 27 H6
q
EOF

#######
# ETOH
gmx make_ndx -f ethanol/md_ETOH_1us.tpr -o ethanol/index_rg.ndx << EOF
r 23 - 28
r 143 - 148
r 88 - 93
r 208 - 213
r 74 - 79
r 194 - 199
20 | 21 | 22 | 23 | 24 | 19
name 25 HAUT
del 19 
del 19 
del 19 
del 19 
del 19 
del 19 
r 33 - 38
r 153 - 158
r 65 - 70
r 185 - 190
r 99 - 104
r 219 - 224
20 | 21 | 22 | 23 | 24 | 25
del 20
del 20
del 20
del 20
del 20
del 20
name 20 MILIEU
r 43 - 48
r 163 - 168
r 55 - 60
r 175 - 180
r 110 - 115
r 230 - 235
26 | 21 | 22 | 23 | 24 | 25
del 21
del 21
del 21
del 21
del 21
del 21
name 21 BAS
r 18 - 49
name 22 H1
r 55 - 77 
name 23 H2
r 84-117
name 24 H3
r 138- 169
name 25 H4
r 175- 197
name 26 H5
r 204-237
name 27 H6
q
EOF

#######
# BOTH
gmx make_ndx -f both/md_ETHE_ETOH_1us.tpr -o both/index_rg.ndx << EOF
r 23 - 28
r 143 - 148
r 88 - 93
r 208 - 213
r 74 - 79
r 194 - 199
20 | 21 | 22 | 23 | 24 | 25
name 26 HAUT
del 20 
del 20 
del 20 
del 20 
del 20 
del 20 
r 33 - 38
r 153 - 158
r 65 - 70
r 185 - 190
r 99 - 104
r 219 - 224
26 | 21 | 22 | 23 | 24 | 25
del 21
del 21
del 21
del 21
del 21
del 21
name 21 MILIEU
r 43 - 48
r 163 - 168
r 55 - 60
r 175 - 180
r 110 - 115
r 230 - 235
26 | 27 | 22 | 23 | 24 | 25
del 22
del 22
del 22
del 22
del 22
del 22
name 22 BAS
r 18 - 49
name 23 H1
r 55 - 77 
name 24 H2
r 84-117
name 25 H3
r 138- 169
name 26 H4
r 175- 197
name 27 H5
r 204-237
name 28 H6
q
EOF

echo "18" | gmx gyrate -f alone/md_VvETR2_1us_fit.xtc -s alone/md_VvETR2_1us.gro -n alone/index_rg.ndx -o alone/gyration_haut.xvg
echo "19" | gmx gyrate -f alone/md_VvETR2_1us_fit.xtc -s alone/md_VvETR2_1us.gro -n alone/index_rg.ndx -o alone/gyration_milieu.xvg
echo "20" | gmx gyrate -f alone/md_VvETR2_1us_fit.xtc -s alone/md_VvETR2_1us.gro -n alone/index_rg.ndx -o alone/gyration_bas.xvg

echo "19" | gmx gyrate -f ethylene/md_ETHE_1us_fit.xtc -s ethylene/md_ETHE_1us.gro -n ethylene/index_rg.ndx -o ethylene/gyration_haut.xvg
echo "20" | gmx gyrate -f ethylene/md_ETHE_1us_fit.xtc -s ethylene/md_ETHE_1us.gro -n ethylene/index_rg.ndx -o ethylene/gyration_milieu.xvg
echo "21" | gmx gyrate -f ethylene/md_ETHE_1us_fit.xtc -s ethylene/md_ETHE_1us.gro -n ethylene/index_rg.ndx -o ethylene/gyration_bas.xvg

echo "19" | gmx gyrate -f ethanol/md_ETOH_1us_fit.xtc -s ethanol/md_ETOH_1us.gro -n ethanol/index_rg.ndx -o ethanol/gyration_haut.xvg
echo "20" | gmx gyrate -f ethanol/md_ETOH_1us_fit.xtc -s ethanol/md_ETOH_1us.gro -n ethanol/index_rg.ndx -o ethanol/gyration_milieu.xvg
echo "21" | gmx gyrate -f ethanol/md_ETOH_1us_fit.xtc -s ethanol/md_ETOH_1us.gro -n ethanol/index_rg.ndx -o ethanol/gyration_bas.xvg

echo "20" | gmx gyrate -f both/md_ETHE_ETOH_1us_fit.xtc -s both/md_ETHE_ETOH_1us.gro -n both/index_rg.ndx -o both/gyration_haut.xvg
echo "21" | gmx gyrate -f both/md_ETHE_ETOH_1us_fit.xtc -s both/md_ETHE_ETOH_1us.gro -n both/index_rg.ndx -o both/gyration_milieu.xvg
echo "22" | gmx gyrate -f both/md_ETHE_ETOH_1us_fit.xtc -s both/md_ETHE_ETOH_1us.gro -n both/index_rg.ndx -o both/gyration_bas.xvg

######################
#### Interactions ####
######################

##################
#### Clusters ####
##################

#gmx cluster -s topol.tpr -f traj_fit.xtc -o clusters.xpm -g cluster.log -cl centers.pdb -method gromos -cutoff 0.2
echo "3" "1" | gmx cluster -s alone/md_VvETR2_1us.tpr -f alone/md_VvETR2_1us_fit.xtc -o alone/clusters.xpm -g alone/cluster.log -cl alone/centers.pdb -method gromos -cutoff 0.2 
echo "3" "1" | gmx cluster -s ethylene/md_ETHE_1us.tpr -f ethylene/md_ETHE_1us_fit.xtc -o ethylene/clusters.xpm -g ethylene/cluster.log -cl ethylene/centers.pdb -method gromos -cutoff 0.2 
echo "3" "1" | gmx cluster -s ethanol/md_ETOH_1us.tpr -f ethanol/md_ETOH_1us_fit.xtc -o ethanol/clusters.xpm -g ethanol/cluster.log -cl ethanol/centers.pdb -method gromos -cutoff 0.2
echo "3" "1" | gmx cluster -s both/md_ETHE_ETOH_1us.tpr -f both/md_ETHE_ETOH_1us_fit.xtc -o both/clusters.xpm -g both/cluster.log -cl both/centers.pdb -method gromos -cutoff 0.2




