
####Run AutogeneS, a gene selection tool for Deconvolution experiments
### The optimizer is run on fixed mode to get specific number of genes. 


import os
import pandas as pd
import anndata
import autogenes as ag

def get_project_root(depth):
    script_path = os.path.abspath(__file__)
    project_root = os.path.join(script_path, *[os.path.pardir] * depth)
    return os.path.abspath(project_root)

def load_csv_data(file_path, index_col='Gene'):
    data = pd.read_csv(file_path, sep=",")
    data.columns = data.iloc[0]  # Set the first row as column headers
    data = data[1:]  # Remove the first row after setting it as the header
    data.index.name = index_col  # Ensure the index name is a string
    return data.set_index(data.columns[0])  # Reassign index to ensure it's properly named

def load_csv_data2(file_path):
    data = pd.read_csv(file_path, header=None, index_col=0)
    data.columns = data.iloc[0]  # Set the first row as column headers
    data = data[1:]
    data.index.name = None  # Ensure the index name is None or a string
    return data

if __name__ == "__main__":
    project_root = get_project_root(depth=4)

    spatial_data = load_csv_data(
        os.path.join(project_root, "Processing/BLADE/Data/Processed_BLADE/BLADE_spatialfinal.csv"),
        index_col='Gene'
    )

    scrna_mean = load_csv_data(
        os.path.join(project_root, "Processing/BLADE/Data/Processed_BLADE/BLADE_scrna_Mean.csv"),
        index_col='Gene'
    )

    scrna_sd = load_csv_data(
        os.path.join(project_root, "Processing/BLADE/Data/Processed_BLADE/BLADE_scrna_SD.csv"),
        index_col='Gene'
    )

    scrna = load_csv_data2(
        os.path.join(project_root, "Processing/BLADE/Data/Processed_BLADE/BLADE_scrna_fullvariable.csv")
    ).T.apply(pd.to_numeric, errors='coerce')


    scrna_meta = load_csv_data2(
        os.path.join(project_root, "Processing/BLADE/Data/Processed_BLADE/BLADE_scrna_meta.csv")
    )

    #variable_features = ["Ly6g6e","Vip","Cck","Trh","S100a5","Gal","Ptgds","Cdhr1","Tac2","Slc17a7","Fst","Igfbp2","Nmb","Apod","Eomes","Gng13","Npy","Th","Tmsb10","Calca","Gm2694","Reln","Ptprd","Olfm1","Id2","Spp1","Kcna2","Ptn","Rab3b","Calb1","Pak1","Tac1","Neurod1","Calb2","Plp1","Arap2","Cacnb4","Ms4a15","Vsnl1","Sv2b","Sncb","Hbb-bs","Plcxd2","Fabp7","Parm1","Pcp4","Gm525","Gabra1","Olfr1295","Pcp2","Sst","Lrrc55","Hba-a1","Hba-a2","Elavl2","Scube1","Apoe","Kitl","C1qa","Shisa3","Tnni3","Samsn1","Ndnf","Gm12027","Lhx1os","Ebf1","Rprm","Ntng1","A230065H16Rik","Crh","Ebf3","Shisa6","Cartpt","Hbb-bt","Yjefn3","Pip5k1b","Olfr1420","Pkib","Syt7","Cbln4","Syndig1l","Mitf","Map1b","C1qc","Nrn1","C1qb","Ankrd34c","Casp1","Coch","Thsd7a","Fxyd7","Cd74","Elavl4","Slc24a2","Tnni1","Apoc1","Dnajc6","Coro6","Zfp804a","Cltc","Car2","Ppm1j","Cbln2","Fxyd1","Hexb","Stmn1","Nefh","Igf1","Cldn5","Ttll7","Rasd2","Gng8","Lyz2","4930523C07Rik","Tbr1","Slc17a6","Edn1","Apba1","Foxp2","Cadps2","Slc20a1","Nppa","Nrsn1","Cpe","Tshz2","Doc2g","Tpbg","6330403K07Rik","Fos","Omp","Ogfrl1","Unc80","Actb","Mgp","Npb","Lhfp","Gap43","Gfy","Adk","Fstl5","Crip1","Ehd3","Rims1","Nrp2","Gucy1a3","AI593442","Xist","Gabrb1","Ywhag","Deptor","Nhlh2","Nxph4","Scn2b","Igfbp5","Hspa4l","Tmem132d","Hlf","Stoml3","Nts","Stmn2","Ebf2","Meg3","Cd24a","Pam","Snrk","Tmem163","Ttpal","Pclo","Nell1","Vtn","Cst3","Nefm","Tpbgl","Cbln1","Fcrls","Prss12","Fkbp1b","Scgn","Crabp1","Grem1","Uchl1","Lor","Lbh","Cux2","Npr1","Nsf","Ctss","Clstn2","Clmp","Lmcd1","Apold1","Zwint","Mpped1","Sytl1","Tnnt1","Mal","Rab37","Cox6a2","Resp18","Sln","Kcng1","Lgals1","Ctgf","Tspan9","Fam171b","Slco1a4","P2ry12","Pfn2","Chgb","Adcyap1","Gpr101","Scn9a","Cpne6","Sult1d1","Nrip3","Serpine2","Pvalb","Nptx1","Rxfp1","Kcnip2","Cldn11","Mdga1","Crim1","Arhgdig","Cplx2","Rassf4","Snhg11","Slc6a11","Tgfa","Hpcal1","Alpk3","Slc6a7","Socs3","Tyrobp","Ermn","Tspan13","Ptprt","Snca","Nr2f1","Higd1b","Snap25","Mapk10","Gad2","Ush1g","Tppp","Hmcn1","Gm14029","Tyrp1","Kcnc1","Ttll6","Olfm3","Rab26","Grm1","Rcan2","Tubb4a","Mslnl","Olfr1288","Slc38a1","Nptx2","Erc2","Ajap1","Sall3","Sfrp5","Lhfpl3","Samd3","Igfbp7","Akap13","Raver2","Chrna4","S100a10","Tnfrsf12a","Atp6ap1l","Rasgrp1","Fam122b","Syt1","Spock1","Ttr","Ralyl","Syp","Nnat","Nt5dc3","Iqgap2","Trp73","Ccl4","Rnf152","Shisa9","Ypel2","Ifitm3","Plxna2","Ccl7","Atp6v1a","Zbtb20","Fam63b","Serpini1","Dlk1","Chrna2","Mgst1","Fat3","Gm1698","Gm27199","Sdc3","Ryr1","Slc6a12","BC048546","Kctd12","Cntnap2","Bhmt","Id3","Rgs5","Crhbp","Tmem40","Pantr1","Ywhah","Cacng7","Atp6v1c2","Ndst4","Tmem37","Prss56","Pgam2","Atp1a2","Tspan2","Kif1b","Cplx1","RP23-407N2.2","Dner","Ghrh","Kit","Nav1","Ndrg4","Osbpl1a","Rab6a","Tnfaip8","Txnrd3","A930019D19Rik","Neurod2","Myh14","Csf1r","Rit2","Car10","3632451O06Rik","Bdnf","Zdbf2","Edil3","Cox8b","Ly6c2","Cryab","Ccl8","Fcer1g","Cdr2","Scn1b","Tfap2c","Slc5a5","Timp2","Pax6","Vwa3b","Mustn1","Gpnmb","Slc35g1","Ccnd2","Kif5c","Strbp","Ntrk3","Etnk1","Inadl","Cst12","Serpinf1","Pou3f1","Jun","Nefl","Sp110","Klf2","Dst","Impact","Nap1l5","F630028O10Rik","Slc47a1","Olfr1297","Cpm","Nnmt","Cx3cr1","Slc9a4","Scn8a","Unc5d","Gja1","Fxyd6","Nr6a1os","Tnfrsf18","Napsa","Cort","Tsku","Wipi1","Dkkl1","Nccrp1","Rbms3","Mpzl2","Chn1","Klc1","Acbd7","Gm38112","Cd9","Bid","Mfge8","Ttbk2","Rgs1","Kcnip4","Lingo1","Lcat","Frzb","Chst15","Nr4a1","Nrgn","Slc38a3","Gng3","Cyfip2","Sorcs3","Dock9","Igfbpl1","Adm","Trf","Nr2f2","Necab2","Ramp3","Olfr1308","Nwd2","Wnt7b","Nrp1","Dmrta2os","Ppef2","Alcam","Tbx21","Gpr26","Sash1","Lynx1","Efhd2","Olfr309","Slc36a2","Mif","Syt4","Stxbp1","Lhx1","Prl","Thy1","Man1a2","Tcap","Chga","C1ql1","Stom","Gm43603","Dsp","Gm27217","Sox11","Prune2","Elovl6","Ctsb","Cdkn1c","Lin7a","Ppfibp1","Tc2n","Ywhab","Faim2","Tpt1","Htr3a","Clec4a3","Tmie","Rgs11","Pcp4l1","Ctxn3","Junb","Kcnc3","Trnp1","Prokr2","Barhl2","Rgcc","Cobl","Laptm5","Mt1","Csmd2","Qrfpr","Tmem132c","Samd5","Rbms1","9330158H04Rik","Fibcd1","Gm9844","Bsn","Vgf","Klf7","Scd1","Gm7233","Ccl9","1700084C06Rik","Dcdc2a","Vldlr","Plac8","Alox8","Adrb2","Zic1","Cnp","Scg2","Hapln4","Fcgr4","Olfr1301","Ms4a4b","Jam2","Pnoc","Arrb1","Chrdl1","Chst8","Rnf14","Brca1","Gpx3","Mllt11","Nlgn2","Hmgcs2","Ppp1r1b","1700019D03Rik","Hapln1","Ptprc","Maf","Sept3","H2-Ab1","1500009C09Rik","Ptchd4","Caln1","Filip1","Ifit2","1110002O04Rik","Olfr1260","Ppp3r2","Pkp3","Sct","Cd3g","Scn4b","Sun3","Spdl1","Fam110c","Olfr209","Slc22a2","H2-Ob","Adgre4","Tmc1","Olfr788","Gm37240","Padi3","RP24-173O20.2","Olfr118","A630095E13Rik","Trem2","Mag","Sparc","Dbpht2","Cyp4x1","Atp5g1","Lypd1","Gls","1500015O10Rik","Cdh4","Sh3gl2","Ier2","Cpne7","Dpysl3","Lmtk2","Ly6a","Itm2a","Slc38a11","Cytip","Rxfp3","Plpp3","Prok2","Adamts1","Dhcr24","Sema3c","Adora1","Slc38a2","Cdh13","Rims3","Cyr61","Rhcg","Mgat4c","Scn2a1","Brinp3","Cd44","5930412G12Rik","Slc5a8","Gm28050","Tmem178b","Hspa12a","Cdkn1b","Itgb2","Pdlim3","Cadm1","Rspo1","Nrn1l","Notch2","Ap3s1","Amhr2","Chchd10","Atrnl1","Ifi27","Camkk2","Vamp2","Tubb3","Fam196a","Cntnap5a","Slc1a3","Lpgat1","Ywhaz","Crip2","Fut9","Nacc2","Inhba","Sepp1","Lmo4","Sema3d","Mgat5","Ifi27l2a","Ttyh2","Gdf10","Lhx9","E530001K10Rik","Fetub","Gm20632","Rbm47","Miat","Sstr2","Trdn","P2ry1","Asap1","C130026I21Rik","Pls1","Pla2g7","Pou3f2","Grm7","Acan","Lsamp","Nfia","St8sia1","Cdkn1a","Ddc","Gm26522","Adap2os","4931423N10Rik","Ly6c1","Fam126b","Hes5","Fcgr3","Pf4","Rell2","Akap7","Hcn1","Syt17","Dynll2","Evi5","Gm11992","Cnga2","Lrrtm1","Selplg","Cd164","AF529169","Fam84b","Epcam","Gpr37l1","Pitpnc1","Hsp90aa1","Snrpn","Scrg1","L1cam","Vdr","Rgs4","Amz1","St8sia6","Kifc3","Islr","Tmem108","Fam183b","Stxbp6","Ptpre","Rasa4","Eps15","Timp3","Ppm1h","Dnm3","Rgs10","Atf5","Cnr1","Slc22a23","Scara5","Slc6a1","Adam11","Oxr1","Otos","Lad1","2600014E21Rik","Ercc6l","Gm42949","Ninj2","C1rl","RP24-118K11.4","Gm9930","Mcemp1","Msr1","A530076I17Rik","Ikzf3","Gsdma","Pfn3","9030624G23Rik","Ly6d","Prss41","Umodl1","Ptprcap","Ttll9","Tmem26","Gm16552","Napb","Lhx2","Hapln3","Scyl2","Kif26b","Bcan","Lst1","Dnaja4","Car4","Siglech","Diras1","A830018L16Rik","Sema3f","Mtor","Mrc1","Rspo3","Pcdh20","Kcnk15","Zfp869","Fam83h","2010005H15Rik","Stmn4","Xrcc3","Nxph1","Dynll1","Tex40","Tnfaip6","Prss23","Olfml3","Nbeal2","Rab15","Adamts7","Ghr","Cdhr3","Npy1r","Reep5","Dnm1","Cysltr1","Gpr183","6530403H02Rik","Cxcl1","RP23-172P1.4","Ido1","Cdcp1","Gpr141","9330161L09Rik","Hs3st6","Aldh1a7","Vsx1","Gm29292","Hgfac","Mypn","Vegfa","Myof","Tbc1d4","Coro1a","Kcnab2","Fam92b","Kcnq3","Hcrtr2","Olfr491","Ptafr","Tenm1","Cebpb","Spock2","Fam168a","Fabp3","Htr2c","Cox7a1","6430548M08Rik","Tcf15","Rac3","Plxnc1","Tuba4a","Htra1","Gsn","Egr1","Abi3bp","Cend1","Gpr34","Nanos2","Wbscr17","Gm11549","Erbb3","Spon1","Rap1gap2","Adgra1","Sema6a","Spata6","Sema7a","Pnldc1","Hopx","Sncg","H2-Eb1","Rtn4rl1","Kif5b","Klhl10","Cxadr","Cxcl14","Srgn","Ptpn4","Mbp","Hs3st4","Nradd","Brsk2","Pdgfrb","Sp8","Mylk","Matn4","Msi2","Creb5","Mrap2","Ppp1r14a","Rmst","Diablo","Tmod2","Pde1c","Gnal","Kcnc2","Sema3e","Ar","Efnb2","Cited4","Mcf2l","Rreb1","Cnrip1","Il1b","Gm4779","A730046J19Rik","Gm128","Cxcl9","8430408G22Rik","Hbb-y","Otx2","Slc17a8","Olfr727","Tmem8c","Olfr1229","Kcns1","Fate1","1700064H15Rik","Gm17660","Bcl3","Zmynd15","AU022754","Sprr1a","1700028P14Rik","Pcdh8","Pcbd1","Rasgrf2","Cngb1","Trpc3","1700012B09Rik","G0s2","Map1lc3a","Tcerg1l","Akap17b","Ly86","Mecom","Itga1","Rgs16","Nupr1","Ramp1","Dclk1","Zfp462","Arpc5l","Smad7","Bcl11b","Klf6","Vwc2l","Spred2","Spidr","Hspb1","Slc7a10","Tbc1d30","Rgs7bp","Cyth1","Gm9925","Igfbp4","Clstn1","Agt","Slco3a1","Rgs8","Stard5","Sema6b","Slc25a22","Whrn","Trpc7","Gm10382","Enc1","Nppc","Gnb1","Ncdn","Flt1","Bhlhe22","Gm20515","Syt11","Golga7b","Prkg2","Mn1","Rasgef1b","Tspan17","Arhgap26","Cpne4","Raph1","Rgs2","Aspa","Pex5l","Emx1","Clstn3","Penk","F3","Nr4a2","Otop2","Cxcr4","Mia","Dpysl2","Tubb4b","Myo5a","1110017D15Rik","Lamp5","Zfp36","Oprk1","Snap91","Nos1","Plk2","Kirrel3","Ngfr","Ap1s3","Fam134b","Gpr88","Rhobtb2","B3galt2","Lgi2","Lypla1","Usp32","Cxcl12","S100a6","Syn1","Ttc39b","Chl1","Erdr1","Slc30a3","Slc6a3","Clec12a","Duox2","Tacr3","2410021H03Rik","Map2k4","Uncx","Clca3a1","Angptl1","Dmrtb1","Zcchc12","Cmtm5","Thsd7b","Gfra2","Nkx6-2","Ephb1","Fmo1","Dbi","Wdr7","Cntn6","Gpr22","Unc13c","Rasl11b","Cycs","1110008P14Rik","Kpnb1","Diras2","Cpeb1","Xk","Cpb1","Tas1r3","Gimap9","Ghrl","Cd163","Styk1","Art4","Fxyd3","G430095P16Rik","Plcg2","Gm3411","Sgcg","S100z","Gm17382","4931440J10Rik","Gm10399","Gpsm1","Bmp2","Nell2","Atp6v1g2","Dusp3","Sfrp2","Marcksl1","Nfib","Map4","Pitpna","Hpgd","Ctsd","Spry2","Dzank1","Dkk3","Clic1","Tmem151a","Rtn3","Rabgap1l","B2m","Pcdh15","Gm3764","Nudt4","Nrtn","Gng2","Ccl3","Trp53i11","Scrn1","Gm42722","Dcun1d4","Kctd2","Gabra3","Camk2n1","Pm20d2","Ano1","Prr15","Nlk","Gm26782","Chrm2","Fryl","Plekhb2","Rassf3","Khdrbs2","Ccdc80","Tmem132a","Zfp503","Cbx3","Cyba","Pltp","AI504432","Slc15a3","Aif1","Snx10","Pja2","Cd302","Kif21a","Igf2","Cntn5","Bhlhe41","Cwh43","Endou","Vash2","Enho","Gpc3","Ralgapa2","Nmnat2","Ugt8a","Kctd14","Atp9a","Iqsec1","Nrxn3","2810030D12Rik","Atp5g3","Lgr5","Gad1","AI662270","Sema3a","Slc12a5","Tiam2","Lcp1","Pth2","Tstd1","Elmod1","Rerg","Elovl4","Il13ra2","Tspan4","Pbk","Tlk1","Pak3","Kcnh8","Pvrl3","Anln","Aldoa","Plxdc1","Tmem266","Ptchd2","Dynlt3","Zfp365","AI413582","Efcab10","Cyp26b1","Gm26759","Pllp","Ppp1r17","Nat8l","Bcl2","Carns1","Atp2b4","Tmem119","Sema5b","Bace2","Mog","Cadps","Pnma2","Sgk1","Fam73b","Itpr1","Sod3","Sorl1","Vat1","Arhgef9","Rims4","Dnah7c","Gm20707","Wee2","Galnt6","Gm16054","Gm26768","E030044B06Rik","Spc25","Casq2","Gprin3","RP23-36D15.6","Ckap2","Il17rd","Itgb7","Cnga4","Mafb","Lgi3","Gpr12","Unc5c","Cerk","Btbd3","Sertm1","Arpp21","Fam73a","Hs3st5","Oxtr","Luzp1","Smim13","Tceb2","Adam33","Nfix","Tram2","Slitrk6","Rasd1","Prkar1b","Rnd1","Unc5b","Gda","Syt6","Mro","Islr2","Selm","Cecr6","Ctxn2","Irf8","Ptprz1","Mrgpre","Pgam1","Sox4","Ly6e","Gabra5","Tubb5","Ctla2a","Lancl3","Prkce","Gbp2","Ubtd2","Uhmk1","Cpne5","Fundc1","Gm3696","Ncald","Garnl3","Pea15a","Cyp27a1","Galnt12","Cyp2a5","2010007H06Rik","Wnt6","Gm26771","Mctp2","Tmsb4x","Mtpn","Fam124a","Abcd3","Clvs2","Dtx1","Col11a1","Peg3","Syt9","Ap3b2","Scoc","Tspyl4","Cldn10","Cdkl2","Aak1","Pid1","Mobp","Vsir","Ifit3","Calm2","Pacsin3","Mdh1","Net1","Disp2","Npepps","Ldhb","Arhgap24","Plcb4","Prkar1a","R3hdm1","Heg1","A930011G23Rik","Gpr37","Serpine3","Fgf1","Cyp39a1","Arhgef4","Map2k1","Pou3f3","G3bp2","Cd48","Sowahb","Hrct1","Emp1","Abca8a","Zbtb4","Fgf12","Tmem74bos","Arid5b","Cfap69","Gm13629","Sidt1","Prkacb","Meig1","Sh3rf1","Btbd17","Ghitm","Ptgs1","Gjd2","Ctnnd2","Clec2g","Bend6","Traf3ip3","Cfap77","Bnc2","Gm12999","Tbx3os1","Gm44029","Kbtbd12","Chek1","Krt20","Mndal","Ace2","Ptgfr","Gm16230","Tex29","Il17b","Ms4a6b","Chia1","Gm43591","Gm29539","Clcnkb","Dyrk4","Upp1","Aspn","Ahnak2","Ptpru","1700097N02Rik","Medag","Nexn","Spef2","Mapk8","Rab11fip5","Cacnb2","Mef2a","Olig1","Eid1","2610528A11Rik","Adarb1","Riiad1","Cib2","Apbb2","Hfe","Gchfr","Rell1","Slco1c1","Grin2c","Ets2","Clasp2","Ifitm10","2310040G07Rik","Zim1","Hcst","Gm26772","Rap1gds1","Ngf","H3f3b","Plekha6","Sfrp1","Cplx3","Rtn4r","Cntnap1","Far2","Psat1","Tmprss5","Eno2","Il33","Neat1","Rad51ap1","Lmo2","Eya2","Cited2","Tmem176b","Adamts2","Irf2bp2","Nlgn1","Rnf11","Rgs7","Trak2","Sparcl1","Spred1","Scn1a","Tmem100","Hmgb2","Prr13","Usp31","Nrarp","Rab3a","Camk4","Pxmp2","Sh3bgrl3","Eif5a2","Rora","Ubash3b","Tubb2b","Pip4k2c","Vim","Gm8104","Myc","Rad54b","Lockd","Plxna4","Syn2","Cdk5r1","Arl6ip1","Srxn1","Cntn4","Abhd3","Celf6","Hist1h4d","Pla2g16","Lmo7","Slc9a3r1","Syt10","Ccdc146","Ppp2r2b","Prnp","Btg2","Rgs3","Acvr1c","6030419C18Rik","Kcnk2","Nme1","D430036J16Rik","Pdcd4","Tmem164","Dnajc5","Wif1","Hmbox1","Gas7","Col23a1","Sdc4","Zic5","RP23-291B1.2","Calm1","Agap1","Mlc1","Tspan15","Hspb3","Cd59a","Abhd11os","Lipa","Cab39","Gabrg2","Foxq1","Trpc5","Dok5","Apba2","Il15ra","Lrrc17","Skap1","Gm7120","4933424G06Rik","Tnfrsf11a","Nusap1","Cybb","Poln","Oas1c","Cftr","Btbd16","Phf11d","Flt4","Ccr10","Itgb3","Gm28557","Zpld1","Hsf2bp","Vwa7","Pirt","Hpse","P2ry6","Cyp1b1","St6galnac5","Nampt","Mapt","Htr5a","Plek","Eif4a2","Gng11","Aldh1a3","Nhs","Mapk3","Fam107b","Lgmn","Insm1","Bcar1","Ppp3ca","Gm17634","Shcbp1l","Agmat","Htra3","Pi16","Piezo2","Retn","Tnfrsf11b","Fbln2","Adcy8","Crhr1","Ddah2","S100a16","Narf","Pvrl1","Wdfy3","Cacna2d3","Ptprk","Pcdh17","Gsk3b","Atf3","Abtb2","Paqr6","1110002E22Rik","Rab6b","App","Wnt10a","Tmem65","Rap1gap","Ppp2r2c","Gm1604a","Pcdh7","2900011O08Rik","Cybrd1","Actg1","Aqp4","Wfs1","Epb41l1","Fbln5","S100b","Grb2","Plcxd3","Limch1","Dlgap1","Slc40a1","Tagln3","Lrp1b","Ifitm2","Col24a1","Klhdc8a","Ass1","Ptx3","Drd1","Necab1","Klf3","Inpp4a","Kcnmb2","Wwc1","Pgm2l1","Dnah7b","Hist2h2be","Brinp1","Rorb","Clu","AI464131","Plin2","Fitm2","Fxyd2","Phgdh","Eif2s3x","Egf","Gatm","Gm16196","Hectd2","Clgn","Adra1b","Robo2","Inpp5a","Grin3a","Celf4","Dcn","Rhoj","Frmpd3","Gramd1b","Homer2","Mapre2","Dio3","Gla","Tnc","AW551984","Tpi1","Tubb6","Glra2","Mtfr1","Aqp1","Plk5","Sept8","AC168977.1","Got2","Basp1","Smap1","Lag3","Plpp4","Prdx5","Hrk","Tsix","Cited1","Galnt14","Hgf","Nkain3","Fam150b","Xlr3a","Aebp1","Rapgef3","Bean1","Arhgap32","Tmem156","Sox6os","2210011C24Rik","Ptpn20","Tekt4","Tlr7","Slc25a24","AI838599","4930515B02Rik","Cxcl10","Nfe2l3","Il12rb2","Gm15704","Klk6","Pidd1","4933406P04Rik","Cd209g","Enpp6","Slc38a8","Gm3667","Arl11","Cd109","Gm12031","Gm15893","Wfikkn2","Smoc2","Guca1b","Cfap53","Nppb","Amph","Lama2","Sh3rf3","Msn","Map7d2","Rfk","Efna5","Kcna6","Tgfb2","Ctxn1","B230312C02Rik","Mbnl1","Pcdh19","Cep170b","Chchd2","Atp6v0a1","Adamtsl1","Bgn","Cfh","Sntb2","Erbb4","Stac","Palmd","H2-Aa","Mfsd6","Galntl6","Ccdc85a","Rpp25","Cnnm1","Sh3bgrl2","Cdk6","Hspb8","Plce1","Slc29a4","4930426D05Rik","Atp6v1b2","Gm26737","Col25a1","2810414N06Rik","Rgs6","Spint2","Ndufa4l2","Enpp2","Pafah1b1","Nrip1","Rgs17","Zfp385b","Pde10a","Sdc2","Lrrc4b","Mrfap1","H2-K1","Lrfn5","Sdcbp","Isg15","Aifm3","Atp5a1","Gstm7","Arl8a","Smad1","Cpd","Rbm34","Ube2q1","Chil1","9130008F23Rik","Pthlh","Smim22","Tgoln1","Ptgs2","Tmem204","Micu3","Sdr39u1","Sphkap","Kcnj12","Rbp4","Postn","Casz1","Slc7a7","Tigar","Herc3","Vdac1","Fam19a2","Fmnl2","2900055J20Rik","Sox10","Prdx6","Hes1","Gm15867","Txnip","Aga","Ppp1r2","Trhr","Fam46a","Pianp","Spdya","Klhl30","6430710C18Rik","Stk33","Vit","Gm13499","Rab25","Tagln","Kdm6b","Chd3os","Cabp1","Strip2","Ahnak","Hspa8","Bex2","Akap12","1700003F12Rik","Ptgis","Gm30173","Ccdc24","Gm15747","Gm26847","Tepp","Adcy4","Mcidas","Gm9920","P2rx6","Upk1b","Ly6g6d","Chp1","Acot7","Tnpo2","Ntrk2","Pkia","Smyd2","Bcat1","Nkd1","Fosl2","Efemp1","Zyg11b","Slc7a1","Rab3gap2","Gpr174","Plekhd1","Cthrc1","Mpeg1","Arc","Sox5","Pde1a","Tmem135","Xylt1","Gpr45","Tm4sf1","Hmgn5","Ngfrap1","Amn","Gng4","Ubl4a","Macrod2","Chmp2b","Dusp8","Egfl7","Myh8","Ppp1r36","Gdf11","Grasp","Elovl7","Cd27","Pde5a","Slc1a2","Morf4l1","Cox6b2","Eif4g2","Tmem255a","Gnl3l","Tnr","Dync1i1","Etv1","Exph5","Chrnb4","Slc6a15","Paip1","Mark1","Vstm2b","Isoc1","A330050F15Rik","Efr3b","Sv2c","Phactr2","Ndrg3","Egflam","Clip3","Grik4","Cd28","Opalin","Mcu","Art3","Rab11a","Dusp14","Egr4","Mtus1","Adam23","Avpi1","C230035I16Rik","Hlcs","Adgrl4","Gm16702","Zfp536","Fat1","Coprs","Atp8b1","Bcl6","Tub","Rhoc","Cacna2d2","Srgap1","Epb41l4a","Plin3","Tmem200a","Rbfox1","Lonrf2","Hey2","Mfap2","Hmgcs1","Mgll","Syngr1","Ccnd1","March1","Opn3","Nid2","Pdia5","Mmp14","Cacng4","Abhd2","Nrg1","Gas2l3","Irs4","Pcdh10","Stau2","Cnot3","Slc30a9","Sacs","Dleu7","Pou2f2","Foxp4","Tshz3","Got1","Lrrk2","Spryd7","Dbp","St18","Rrm2","Lum","Lrig3","Ajuba","Dock5","Man2c1os","Gm12279","Abca6","F12","Gpr68","Tmem30a","Rheb","Shh","Rnd3","Osbpl6","C1ql3","Adamts5","Prkar2b","Dock1","Srpk2","Hspa2","Cdc42ep3","Gpr6","2610028E06Rik","Serpine1","Plekhg2","Dqx1","Bfsp2","Fam111a","Grm2","Rec8","Col6a1","Lingo2","Cltb","Rab3c","Btg1","A830010M20Rik","Evi2a","Wnt7a","Eva1c","Mlf1","Pcdh11x","Zfp618","Cyp11a1","Sdk2","Slitrk2","5730403I07Rik","Vwa5b2","Fam101b","Fam159b","Esm1","Arl4c","Kctd8","Lig4","Glrb","Dync1li1","Pglyrp1","Nabp1","Mbnl2","Stx12","Fads6","Mt2","Kifap3","Sox9","Zic4","Fam131c","Glra1","Cacng6","Eef1a2","Ctnnbip1","Gabarapl1","Sema6d","Cpa6","1700016K19Rik","Sema6c","1700086L19Rik","Fam19a1","Lrriq1","Vrk2","Adap2","H2-D1","Sacm1l","Cox5a","Hspa4","Tuba8","Nt5dc2","Lrrc49","Smpdl3a","Igf1r","Fam19a3","Synpo2","Llgl2","Galnt13","Ccng1","Tmem88b","Gm6297","Golim4","Insig1","Fam184a","Npnt","Mapre1","Krt222","Ncan","Crtac1","Nfatc1","9430020K01Rik","Zfp36l1","Pnpo","C330027C09Rik","Rab33a","Kif2a","Optc","Ifi203","Fap","Hsd3b3","Grhl3","Draxin","BB365896","Gm26789","Adamtsl5","1700016D06Rik","Agrp","B230110C06Rik","1810041H14Rik","Tbx18","Gm28111","2310075C17Rik","Slfn2","Npw","Pcdhga4","Ch25h","Filip1l","Ppp1r14c","Mtss1","Amigo2","Card10","Ube2m","Ap1s2","Cd68","Mettl5os","Larp1","Cntn2","Hras","Rragc","Anxa5","Adgrf5","Ppp2ca","Emcn","Pcolce2","Ephx4","A330102I10Rik","Slc7a11","Dusp2","Ripk2","3110082J24Rik","Slc6a17","Fabp5","Runx1t1","Ppp2cb","Gpc6","1810011O10Rik","Rnf17","Thbs1","Ccno","Npr3","Npas4","Limd2","2310022B05Rik"]





    # Create an AnnData object
    adata = anndata.AnnData(scrna)

    # Assign metadata
    adata.obs = scrna_meta

    # Name the index of .var to 'gene_ids'
    adata.var.index.name = 'gene_ids'
    adata.var_names_make_unique()

    # Now safe to call this method
    #adata = adata[:, variable_features]

    print(adata.shape, "adata object shape")

    # Initialize AutoGeneS with the AnnData object
    ag.init(adata, celltype_key='celltypes')  # Adjust 'celltype_key' as necessary

    # Run the optimizer
    ag.optimize(ngen=5000,offspring_size =100, seed=0,nfeatures= 600 , mode = 'fixed')


    # Select the first solution

    #ag.plot(weights=(-1, 0))

    # Extract the selected genes

    selected_genes = ag.select(close_to=(1, 110))

    #selected_genes = ag.select()




    print("Selected genes:", selected_genes)

    selected_gene_names = adata.var_names[selected_genes].tolist()
    print("Selected gene names:", selected_gene_names)

    filename = "Results/selected_genes600.txt"

    # Open the file in write mode and save the list
    with open(filename, 'w') as file:
        for gene in selected_gene_names:
            file.write(f"{gene}\n")

    print(f"Selected gene names saved to {filename}")



