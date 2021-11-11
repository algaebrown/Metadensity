#!/bin/bash
#PBS -q home-yeo
#PBS -N ENCODE_hg38_rmats
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00
#PBS -o /home/hsher/cluster_msg/encode_rmats_EVAN.out
#PBS -e /home/hsher/cluster_msg/encode_rmats_EVAN.err
#PBS -t 0-99


module load rmats/3.2.5
outdir=/projects/ps-yeolab5/encode/rnaseq/alt_splicing_hg38/
basedir=/projects/ps-yeolab5/encode/rnaseq/shrna_knockdown/

bamones=(ENCFF134OIU ENCFF398WLJ ENCFF253HMW ENCFF767VDD ENCFF547ZBH ENCFF259ZYJ ENCFF539VZC ENCFF560FIY ENCFF110BVT ENCFF665WMD ENCFF662XNN ENCFF670UJI ENCFF731OEM ENCFF202OKI ENCFF512ELZ ENCFF611HGV ENCFF205KBS ENCFF384PMZ ENCFF328GAZ ENCFF422YHM ENCFF622ULI ENCFF198BST ENCFF606LIW ENCFF217CYB ENCFF184EXQ ENCFF921OOX ENCFF751ITZ ENCFF771AHF ENCFF144WCH ENCFF324GVC ENCFF185GLF ENCFF891JDX ENCFF992TTP ENCFF973BAQ ENCFF842AXS ENCFF491XQX ENCFF350YQV ENCFF843KJF ENCFF472JQQ ENCFF945HAL ENCFF271ZYT ENCFF721DJW ENCFF754GHW ENCFF742IUO ENCFF555VXW ENCFF468TYN ENCFF870KFE ENCFF449XHR ENCFF815NQY ENCFF201FDG ENCFF631HEW ENCFF997NWS ENCFF652GWU ENCFF079UCH ENCFF927NOX ENCFF842IKX ENCFF650ZEJ ENCFF895YXH ENCFF260OUK ENCFF772JIO ENCFF036KPI ENCFF320LOZ ENCFF126ZCQ ENCFF402XUP ENCFF072KAU ENCFF130KNV ENCFF516KCI ENCFF234KAR ENCFF268SDJ ENCFF081IUU ENCFF055EEK ENCFF357UJY ENCFF458KRT ENCFF410LSF ENCFF244ZQM ENCFF514OEF ENCFF716YLR ENCFF459YCM ENCFF454VQW ENCFF202RTM ENCFF820TIZ ENCFF714UIW ENCFF383DKD ENCFF513STE ENCFF995JDI ENCFF033BYF ENCFF652GQO ENCFF996LET ENCFF519AHH ENCFF885WBL ENCFF928GEG ENCFF425BAT ENCFF467QHA ENCFF314LOA ENCFF324AFT ENCFF995EQN ENCFF985QNQ ENCFF806UPC ENCFF940WSQ ENCFF984MME)

bamtwos=(ENCFF544VWD ENCFF826NMS ENCFF052UJN ENCFF039EFY ENCFF574CZD ENCFF061UBO ENCFF596NAC ENCFF224TGR ENCFF115THD ENCFF062UFO ENCFF747RMM ENCFF850EIP ENCFF874STK ENCFF822TUL ENCFF695SSC ENCFF233YFQ ENCFF614HPB ENCFF642PWD ENCFF275GXI ENCFF471PNG ENCFF472EMZ ENCFF124CQD ENCFF058MXX ENCFF625CKI ENCFF745OSA ENCFF798RSV ENCFF912FLN ENCFF823IDH ENCFF756HIE ENCFF743XCK ENCFF899HDD ENCFF181APT ENCFF072KHY ENCFF991PNX ENCFF543CIU ENCFF998XFR ENCFF432OJL ENCFF624NBF ENCFF663VXW ENCFF730SRO ENCFF574CEC ENCFF345BNB ENCFF338PHG ENCFF206UTF ENCFF533SIN ENCFF593VAX ENCFF755BAZ ENCFF187PGH ENCFF721YQT ENCFF271ABA ENCFF592GCC ENCFF934NGG ENCFF038NBJ ENCFF911NYK ENCFF658DCV ENCFF799XRW ENCFF431LMX ENCFF411QXQ ENCFF122BHD ENCFF940WWC ENCFF964OHS ENCFF781LCV ENCFF641OYB ENCFF580BXZ ENCFF623UWZ ENCFF428OGQ ENCFF752QBE ENCFF285ELE ENCFF208TBN ENCFF319MJT ENCFF838XHX ENCFF039VBF ENCFF569QEZ ENCFF722ZRT ENCFF889FRC ENCFF771HJQ ENCFF315XSA ENCFF541IMJ ENCFF761YCL ENCFF624WAO ENCFF894TVP ENCFF562IQV ENCFF466CKO ENCFF369SFS ENCFF068LRE ENCFF167OCW ENCFF421LNM ENCFF745XZP ENCFF500CQP ENCFF998QZZ ENCFF629DXH ENCFF560KLF ENCFF051CAG ENCFF837FNR ENCFF070OJW ENCFF647IUR ENCFF968NXH ENCFF396WPA ENCFF860DAE ENCFF740UYB)

ctrlones=(ENCFF185IAD ENCFF124BQH ENCFF994YDH ENCFF839SZZ ENCFF310OHO ENCFF694JWV ENCFF371TBZ ENCFF839SZZ ENCFF488WGA ENCFF396YDD ENCFF456VWU ENCFF958NOA ENCFF694JWV ENCFF078LXU ENCFF289TTH ENCFF568DGW ENCFF839SZZ ENCFF549AUI ENCFF925CTT ENCFF529TIM ENCFF951ZGS ENCFF420VMT ENCFF568DGW ENCFF519WJE ENCFF420VMT ENCFF958NOA ENCFF925CTT ENCFF420VMT ENCFF958NOA ENCFF549AUI ENCFF636FCJ ENCFF155ZJN ENCFF994YDH ENCFF549AUI ENCFF665OLB ENCFF232HCC ENCFF568DGW ENCFF176FHW ENCFF396YDD ENCFF376WQN ENCFF549AUI ENCFF556NPQ ENCFF456VWU ENCFF529TIM ENCFF074QTM ENCFF925CTT ENCFF232HCC ENCFF051CXY ENCFF925CTT ENCFF694JWV ENCFF928FXN ENCFF185IAD ENCFF808WAQ ENCFF185IAD ENCFF839SZZ ENCFF839SZZ ENCFF185IAD ENCFF880DNL ENCFF951ZGS ENCFF456VWU ENCFF176FHW ENCFF095ZAS ENCFF420VMT ENCFF556NPQ ENCFF568DGW ENCFF925CTT ENCFF694JWV ENCFF289WRI ENCFF839SZZ ENCFF212WWN ENCFF396YDD ENCFF232HCC ENCFF839SZZ ENCFF636CZX ENCFF456VWU ENCFF636FCJ ENCFF371TBZ ENCFF583JDB ENCFF420VMT ENCFF694JWV ENCFF086SXC ENCFF310OHO ENCFF568DGW ENCFF556NPQ ENCFF958NOA ENCFF051CXY ENCFF095ZAS ENCFF590QAL ENCFF371TBZ ENCFF808WAQ ENCFF124BQH ENCFF396YDD ENCFF289TTH ENCFF964NLY ENCFF844FNN ENCFF808WAQ ENCFF880DNL ENCFF289TTH ENCFF808WAQ ENCFF095ZAS)

ctrltwos=(ENCFF021NQD ENCFF311TMC ENCFF574REP ENCFF161EYD ENCFF743GZQ ENCFF902YHV ENCFF403KGR ENCFF161EYD ENCFF909UYT ENCFF867NHV ENCFF565WRF ENCFF202ZWS ENCFF902YHV ENCFF669EEP ENCFF284NQQ ENCFF691JCT ENCFF161EYD ENCFF761SAC ENCFF360UZH ENCFF361KUQ ENCFF008VSN ENCFF501COM ENCFF691JCT ENCFF272FNP ENCFF501COM ENCFF202ZWS ENCFF360UZH ENCFF501COM ENCFF202ZWS ENCFF761SAC ENCFF936ZIY ENCFF601IUN ENCFF574REP ENCFF761SAC ENCFF615KJS ENCFF979BGR ENCFF691JCT ENCFF122VBO ENCFF867NHV ENCFF231QDM ENCFF761SAC ENCFF598ZGR ENCFF565WRF ENCFF361KUQ ENCFF046OMI ENCFF360UZH ENCFF979BGR ENCFF775HZV ENCFF360UZH ENCFF902YHV ENCFF279WLO ENCFF021NQD ENCFF181IAL ENCFF021NQD ENCFF161EYD ENCFF161EYD ENCFF021NQD ENCFF617TFJ ENCFF008VSN ENCFF565WRF ENCFF122VBO ENCFF782JJL ENCFF501COM ENCFF598ZGR ENCFF691JCT ENCFF360UZH ENCFF902YHV ENCFF197UJB ENCFF161EYD ENCFF698JGB ENCFF867NHV ENCFF979BGR ENCFF161EYD ENCFF737ISE ENCFF565WRF ENCFF936ZIY ENCFF403KGR ENCFF816LBP ENCFF501COM ENCFF902YHV ENCFF652TWL ENCFF743GZQ ENCFF691JCT ENCFF598ZGR ENCFF202ZWS ENCFF775HZV ENCFF782JJL ENCFF286MAR ENCFF403KGR ENCFF181IAL ENCFF311TMC ENCFF867NHV ENCFF284NQQ ENCFF359FKE ENCFF787CCQ ENCFF181IAL ENCFF617TFJ ENCFF284NQQ ENCFF181IAL ENCFF782JJL)

accessions=(ENCSR896CFV ENCSR560RSZ ENCSR300IEW ENCSR408SDL ENCSR448JAM ENCSR148MQK ENCSR509LIV ENCSR907UTB ENCSR040WAK ENCSR201WFU ENCSR656DQV ENCSR295XKC ENCSR939ZRA ENCSR149DMY ENCSR278CHI ENCSR762FEO ENCSR715XZS ENCSR208GPE ENCSR850CKU ENCSR244SIO ENCSR466NRW ENCSR728BOL ENCSR379VXW ENCSR155EZL ENCSR354XQY ENCSR010ZMZ ENCSR478FJK ENCSR880DEH ENCSR906WTM ENCSR577OVP ENCSR891AXF ENCSR060KRD ENCSR778AJO ENCSR164TLB ENCSR118XYK ENCSR237IWZ ENCSR322XVS ENCSR995RPB ENCSR792XFP ENCSR116QBU ENCSR079LMZ ENCSR546MBH ENCSR047VPW ENCSR967QNT ENCSR454KYR ENCSR016IDR ENCSR116YMU ENCSR047IUS ENCSR313CHR ENCSR017PRS ENCSR518JXY ENCSR527QNC ENCSR529MBZ ENCSR925SYZ ENCSR286OKW ENCSR850PWM ENCSR490DYI ENCSR517JDK ENCSR184YDW ENCSR946OFN ENCSR547NWD ENCSR945GUR ENCSR871BXO ENCSR119QWQ ENCSR598YQX ENCSR710NWE ENCSR597XHH ENCSR385KOY ENCSR783YSQ ENCSR778WPL ENCSR792CBM ENCSR681SMT ENCSR973QSV ENCSR755KOM ENCSR028YAQ ENCSR004RGI ENCSR695XOD ENCSR207QGW ENCSR784FTX ENCSR094KBY ENCSR936TED ENCSR624OUI ENCSR921KDS ENCSR577XBW ENCSR028ITN ENCSR048BWH ENCSR812TLY ENCSR296ERI ENCSR905HID ENCSR902WSK ENCSR117WLY ENCSR545AIK ENCSR840QOH ENCSR562CCA ENCSR599PXD ENCSR077BPR ENCSR958KSY ENCSR409CSO ENCSR222ABK ENCSR494UDF)

lengths=(100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 101 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 101 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 101 100 100 100)
################## AS replicate ####################



bam1=${bamones[$PBS_ARRAYID]}
bam2=${bamtwos[$PBS_ARRAYID]}
ctrl1=${ctrlones[$PBS_ARRAYID]}
ctrl2=${ctrltwos[$PBS_ARRAYID]}
acc=${accessions[$PBS_ARRAYID]}
len=${lengths[$PBS_ARRAYID]}

echo $ctrl2

fullout=$outdir$acc
mkdir $fullout

bamall=$basedir$bam1.bam,$basedir$bam2.bam
ctrlall=$basedir$ctrl1.bam,$basedir$ctrl2.bam

rmats \
-b1 $bamall \
-b2 $ctrlall \
-o $fullout \
-t paired \
-len $len \
-gtf /home/hsher/gencode_coords/gencode.v33.annotation.gtf \
-libType fr-firststrand