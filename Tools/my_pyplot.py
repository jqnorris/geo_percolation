import matplotlib as plt
import matplotlib.pyplot as plt

myRed=(229/255.0,25/255.0,50/255.0)
myBlue=(25/255.0,178/255.0,255/255.0)
myGreen=(50/255.0,255/255.0,0/255.0)
myOrange=(255/255.0,127/255.0,0/255.0)
myViolet=(101/255.0,76/255.0,255/255.0)

plt.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{amsmath}',
       r'\usepackage{textgreek}',
       r'\usepackage{amssymb}',
       r'\usepackage{dejavu}',    # set the normal font here
       r'\renewcommand*\familydefault{\sfdefault}', # Only if the base font of the document is to be typewriter style
       r'\usepackage[LGR, T1]{fontenc}',
       r'\renewcommand{\Gamma}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{71}}}'
       r'\renewcommand{\Delta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{68}}}'
       r'\renewcommand{\Theta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{74}}}'
       r'\renewcommand{\Lambda}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{76}}}'
       r'\renewcommand{\Xi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{88}}}'
       r'\renewcommand{\Pi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{80}}}'
       r'\renewcommand{\Sigma}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{83}}}'
       r'\renewcommand{\Upsilon}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{85}}}'
       r'\renewcommand{\Phi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{70}}}'
       r'\renewcommand{\Psi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{89}}}'
       r'\renewcommand{\Omega}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{87}}}'
       r'\renewcommand{\alpha}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{97}}}'
       r'\renewcommand{\beta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{98}}}',
       r'\renewcommand{\gamma}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{103}}}',
       r'\renewcommand{\delta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{100}}}',
       r'\renewcommand{\epsilon}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{101}}}',
       r'\renewcommand{\zeta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{122}}}',
       r'\renewcommand{\eta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{104}}}',
       r'\renewcommand{\theta}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{106}}}',
       r'\renewcommand{\iota}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{105}}}',
       r'\renewcommand{\kappa}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{107}}}',
       r'\renewcommand{\lambda}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{108}}}',
       r'\renewcommand{\mu}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{109}}}',
       r'\renewcommand{\nu}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{110}}}',
       r'\renewcommand{\xi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{120}}}',
       r'\newcommand{\omicron}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{113}}}',
       r'\renewcommand{\pi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{112}}}',
       r'\renewcommand{\rho}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{114}}}',
       r'\renewcommand{\sigma}{\text{\textsigma}}',
       r'\renewcommand{\tau}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{116}}}',
       r'\renewcommand{\upsilon}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{117}}}',
       r'\renewcommand{\phi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{102}}}',
       r'\renewcommand{\chi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{113}}}',
       r'\renewcommand{\psi}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{121}}}',
       r'\renewcommand{\omega}{\text{\usefont{LGR}{\sfdefault}{m}{it}\symbol{119}}}',
       r'\usepackage[]{sansmath}',  
       r'\sansmath'
] 
