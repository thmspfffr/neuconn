%% nc_allsubj_start.m
% contains all neuconn subject codes

function allsubj = nc_allsubj_start(m)

if m == 1
  
  allsubj{1} = ...
         {'pat01_PPSKIF'     ,1    ,1;...
          'pat02_ONEMIJ'     ,2    ,1;...
          'pat03_TJRSMM'     ,1    ,0;... %
          'pat04_KHIKPR'     ,2    ,0;...
          'pat05_OUIVGU'     ,1    ,1;...
          'pat06_GOWRQM'     ,1    ,1;...
          'pat07_RNJWPJ'     ,2    ,1;...
          'pat08_XXXXXX'     ,2    ,0;...
          'pat09_IPVPHK'     ,1    ,0;...
          'pat10_TEHNFJ'     ,2    ,1;...        
          'pat11_GISNLK'     ,1    ,1;... % SLOW ARTIFACT!
          'pat12_QKPNGN'     ,2    ,1;...
          'pat13_IMHKRQ'     ,2    ,1;...
          'pat14_JRJVQK'     ,2    ,1;...
          'pat15_CFPUMQ'     ,1    ,1;...
          'pat16_JPHDKR'     ,2    ,1;...
          'pat17_OKVSRO'     ,1    ,1;...        
          'pat18_GQPOQJ'     ,2    ,1;...
          'pat19_PMKRRG'     ,1    ,1;...
          'pat20_LNREIL'     ,2    ,1;...
          'pat21_FPPIKN'     ,1    ,1;... % HEAVY MUSCLE ARTIFACTS
          'pat22_QRUECO'     ,2    ,0;... % STRONG EXTERNAL ARTIFACT
          'pat23_NMJQFS'     ,2    ,0;...
          'pat24_HKDILT'     ,2    ,1;...        
          'pat25_MTTKKG'     ,1    ,1;...
          'pat26_XXXXXX'     ,2    ,0;...
          'pat27_XXXXXX'     ,1    ,0;...
          'pat28_QIDKDQ'     ,2    ,1;...
          'pat29_HMJMEK'     ,1    ,1;...
          'pat30_FKHDJP'     ,2    ,1;...
          'pat31_QLJOPI'     ,1    ,1;...
          'pat32_RSPLKE'     ,2    ,1;...
          'pat33_SENQLM'     ,1    ,1;...
          'pat34_XXXXXX'     ,2    ,0;...
          'pat35_EMEMEG'     ,1    ,1;...
          'pat36_ORKLIT'     ,2    ,1;...
          'pat37_FROOIR'     ,2    ,1;...
          'pat38_KQJPMQ'     ,2    ,1;...
          'pat39_VNLSDI'     ,1    ,1;...
          'pat40_FMKLJR'     ,2    ,1;...
          'pat41_SOSNOP'     ,1    ,1;...
          'pat42_KODLSI'     ,2    ,1;...
          'pat43_TRSPIO'     ,1    ,1;...
          'pat44_TTWKMX'     ,1    ,1;...
          'pat45_MLTQLJ'     ,1    ,1;...
          'pat46_JHKSBM'     ,2    ,1;...
          'pat47_IKRBSB'     ,1    ,1;...
          'pat48_PQPNJT'     ,2    ,1;...
          'pat49_LLBPJI'     ,1    ,1};

  allsubj{2} = ...
    {'con01_IIWFGW'    ,1    ,1;...
    'con02_HLJWGC'     ,2    ,1;...
    'con03_XXXXXX'     ,1    ,0;...
    'con04_XXXXXX'     ,2    ,0;...
    'con05_NFLQLC'     ,1    ,1;...
    'con06_QJNTOG'     ,2    ,1;...
    'con07_MPSLOP'     ,1    ,1;...
    'con08_TPKKHG'     ,2    ,1;...
    'con09_XXXXXX'     ,1    ,0;...
    'con10_LEGIKB'     ,2    ,1;...
    'con11_LEBQKB'     ,1    ,1;...
    'con12_IFSBSI'     ,2    ,1;...
    'con13_QJNTOG'     ,2    ,1;...
    'con14_RLODRC'     ,2    ,1;...
    'con15_NJBQSK'     ,1    ,1;... % TASK ORDER???
    'con16_KKJDGL'     ,2    ,1;...
    'con17_QKHTPH'     ,1    ,1;...
    'con18_KDPOIF'     ,2    ,1;...
    'con19_GSJQGL'     ,1    ,1;...
    'con20_RTPSQW'     ,2    ,1;...
    'con21_PLTILJ'     ,1    ,1;...
    'con22_XXXXXX'     ,2    ,0;...
    'con23_VBSOOE'     ,1    ,1;...
    'con24_MKLGNE'     ,2    ,1;...
    'con25_TALVGC'     ,1    ,1;...
    'con26_XXXXXX'     ,2    ,0;...
    'con27_XXXXXX'     ,1    ,0;...
    'con28_GDRXQG'     ,2    ,1;...
    'con29_PMKJIP'     ,1    ,1;...
    'con30_JDQMIL'     ,2    ,1;...
    'con31_PNDMKU'     ,1    ,1;...
    'con32_JITIHV'     ,2    ,1;...
    'con33_JQGIIM'     ,1    ,1;...
    'con34_XXXXXX'     ,2    ,0;...
    'con35_OVLQHN'     ,1    ,1;...
    'con36_NJKVQJ'     ,1    ,1;...
    'con37_LOESLV'     ,1    ,1;...
    'con38_RSJEVI'     ,2    ,1;...
    'con39_GYLBVV'     ,1    ,1;...
    'con40_JKNXJM'     ,2    ,1;...
    'con41_RRPNPO'     ,1    ,1;...
    'con42_NNDLMK'     ,2    ,1;...
    'con43_MNDMLK'     ,1    ,1;...
    'con44_NPDPMM'     ,2    ,1;...
    'con45_GDRXQG'     ,1    ,1;...
    'con46_VVSMVJ'     ,2    ,1;...
    'con47_JOJLGN'     ,1    ,1;...
    'con48_OWJHNM'     ,2    ,1;...
    'con49_SULTRR'     ,1    ,1};

else
  
  allsubj{1} = ...
       {'pat01_PPSKIF'     ,2    ,1;...
        'pat02_ONEMIJ'     ,1    ,1;...
        'pat03_TJRSMM'     ,2    ,0;... %
        'pat04_KHIKPR'     ,1    ,0;...
        'pat05_OUIVGU'     ,2    ,1;...
        'pat06_GOWRQM'     ,2    ,1;...
        'pat07_RNJWPJ'     ,2    ,1;...
        'pat08_XXXXXX'     ,1    ,0;...
        'pat09_IPVPHK'     ,2    ,0;...
        'pat10_TEHNFJ'     ,1    ,1;...        
        'pat11_GISNLK'     ,2    ,1;... % SLOW ARTIFACT!
        'pat12_QKPNGN'     ,1    ,1;...
        'pat13_IMHKRQ'     ,1    ,1;...
        'pat14_JRJVQK'     ,1    ,1;...
        'pat15_CFPUMQ'     ,2    ,1;...
        'pat16_JPHDKR'     ,1    ,1;...
        'pat17_OKVSRO'     ,2    ,1;...        
        'pat18_GQPOQJ'     ,1    ,1;...
        'pat19_PMKRRG'     ,2    ,1;...
        'pat20_LNREIL'     ,1    ,1;...
        'pat21_FPPIKN'     ,2    ,1;... % HEAVY MUSCLE ARTIFACTS
        'pat22_QRUECO'     ,1    ,0;... % STRONG EXTERNAL ARTIFACT
        'pat23_XXXXXX'     ,2    ,0;...
        'pat24_HKDILT'     ,1    ,1;...        
        'pat25_MTTKKG'     ,2    ,1;...
        'pat26_EKKFOI'     ,1    ,1;...
        'pat27_XXXXXX'     ,2    ,0;...
        'pat28_QIDKDQ'     ,1    ,1;...
        'pat29_HMJMEK'     ,1    ,1;...
        'pat30_FKHDJP'     ,1    ,1;...
        'pat31_QLJOPI'     ,2    ,1;...
        'pat32_RSPLKE'     ,1    ,1;...
        'pat33_SENQLM'     ,2    ,1;...
        'pat34_XXXXXX'     ,1    ,0;...
        'pat35_EMEMEG'     ,2    ,1;...
        'pat36_ORKLIT'     ,1    ,1;...
        'pat37_XXXXXX'     ,2    ,0;...
        'pat38_KQJPMQ'     ,1    ,1;...
        'pat39_VNLSDI'     ,2    ,1;...
        'pat40_FMKLJR'     ,1    ,1;...
        'pat41_XXXXXX'     ,2    ,0;...
        'pat42_KODLSI'     ,1    ,1;...
        'pat43_XXXXXX'     ,2    ,0;...
        'pat44_TTWKMX'     ,1    ,1;...
        'pat45_MLTQLJ'     ,2    ,1;...
        'pat46_JHKSBM'     ,1    ,1};

  allsubj{2} = ...
    {'con01_XXXXXX'    ,2    ,1;...
    'con02_HLJWGC'     ,1    ,1;...
    'con03_XXXXXX'     ,2    ,1;...
    'con04_XXXXXX'     ,1    ,1;...
    'con05_NFLQLC'     ,2    ,1;...
    'con06_QJNTOG'     ,1    ,1;...
    'con07_MPSLOP'     ,2    ,1;...
    'con08_TPKKHG'     ,1    ,1;...
    'con09_XXXXXX'     ,2    ,1;...
    'con10_LEGIKB'     ,1    ,1;...
    'con11_LEBQKB'     ,2    ,1;...
    'con12_IFSBSI'     ,1    ,1;...
    'con13_QJNTOG'     ,2    ,1;...
    'con14_RLODRC'     ,1    ,1;...
    'con15_NJBQSK'     ,2    ,1;... % TASK ORDER???
    'con16_XXXXXX'     ,1    ,1;...
    'con17_QKHTPH'     ,2    ,1;...
    'con18_KDPOIF'     ,1    ,1;...
    'con19_KEFELD'     ,2    ,1;...
    'con20_XXXXXX'     ,1    ,1;...
    'con21_PLTILJ'     ,2    ,1;...
    'con22_XXXXXX'     ,1    ,1;...
    'con23_VBSOOE'     ,2    ,1;...
    'con24_MKLGNE'     ,1    ,1;...
    'con25_TALVGC'     ,2    ,1;...
    'con26_XXXXXX'     ,1    ,1;...
    'con27_XXXXXX'     ,2    ,1;...
    'con28_GDRXQG'     ,1    ,1;...
    'con29_PMKJIP'     ,2    ,1;...
    'con30_JDQMIL'     ,1    ,1;...
    'con31_XXXXXX'     ,2    ,1;...
    'con32_JITIHV'     ,1    ,1;...
    'con33_JQGIIM'     ,2    ,1;...
    'con34_XXXXXX'     ,1    ,1;...
    'con35_OVLQHN'     ,2    ,1;...
    'con36_XXXXXX'     ,1    ,1;...
    'con37_XXXXXX'     ,2    ,1;...
    'con38_RSJEVI'     ,1    ,1;...
    'con39_XXXXXX'     ,2    ,1;...
    'con40_JKNXJM'     ,1    ,1};
  
end
