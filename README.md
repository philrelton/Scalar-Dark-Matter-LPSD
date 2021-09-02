# Scalar-Dark-Matter-LPSD
This repository contains code as described in Vermulen et al, 2021. This includes:
- An algorithm for calculating a Logarithmically spaced Power-Spectral-Density (LPSD)
- Examples for how to run the LPSD algorithm
- Code for adding dark matter signals to data for testing. This allows for blind testing

## Requirements:
The file `Scalar-Dark-Matter-LPSD.yml` contains the information required to create a working conda 
environment for this project. The necessary packages can be installed using the command:
`conda env create -f Scalar-Dark-Matter-LPSD.yml`

## Licensing:
The LPSD is a modified version of the version found here https://gitlab.aei.uni-hannover.de/geoq/lpsd.
The original version is unlicensed as such we have requested and recieved written permission from the 
authors to use and adapt their code. We have also recieved permission to license our version under the
GNU GENERAL PUBLIC LICENSE. You can find this license in the LICENSE file in this repository.

If you do publish work using this code, please cite the following papers:
- Tröbs, M. and Heinzel, G., 2006. Improved spectrum estimation from digitized time series on a logarithmic frequency axis. Measurement, 39(2), pp.120-129.
https://www.sciencedirect.com/science/article/pii/S026322410500117X?casa_token=WcRBlyyEABYAAAAA:TqUNIcSN2qWlFMFr1eHROqnafKYUFQq14yDCpJX6S8PE593F9P5LSSOCL4AxL90fxb3PR9gHJw
```
@article{trobs2006improved,
  title={Improved spectrum estimation from digitized time series on a logarithmic frequency axis},
  author={Tr{\"o}bs, Michael and Heinzel, Gerhard},
  journal={Measurement},
  volume={39},
  number={2},
  pages={120--129},
  year={2006},
  publisher={Elsevier}
}
```
- Tröbs, M. and Heinzel, G., 2009. Corrigendum to “Improved spectrum estimation from digitized time series on a logarithmic frequency axis” [Measurement 39 (2006) 120–129]
https://www.sciencedirect.com/science/article/pii/S0263224108000705
```
@article{trobs2009improved,
  title={Improved spectrum estimation from digitized time series on a logarithmic frequency axis (vol 39, pg 120, 2006)},
  author={Tr{\"o}bs, Michael and Heinzel, Gerhard},
  journal={Measurement},
  volume={42},
  number={1},
  pages={170--170},
  year={2009}
}
```
- Vermeulen et al. 2021: https://arxiv.org/pdf/2103.03783.pdf
```
@article{vermeulen2021direct,
  title={Direct limits for scalar field dark matter from a gravitational-wave detector},
  author={Vermeulen, Sander M and Relton, Philip and Grote, Hartmut and Raymond, Vivien and Affeldt, Christoph and Bergamin, Fabio and Bisht, Aparna and Brinkmann, Marc and Danzmann, Karsten and Doravari, Suresh and others},
  journal={arXiv preprint arXiv:2103.03783},
  year={2021}
}
```

Please feel free to use all non-LPSD code here without request. However, if you do, please cite Vermeulen et al.
