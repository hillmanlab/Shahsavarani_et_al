# Shahsavarani_et_al

**GENERAL INFORMATION**  
This data is described in the following publication:  
  
  __Cortex-wide neural dynamics predict behavioral states and provide a neural basis for resting-state dynamic functional connectivity__, Somayeh Shahsavarani1,2,5, David N. Thibodeaux1,5, Weihao Xu1, Sharon H. Kim1, Fatema Lodgher1, Chinwendu Nwokeabia1, Morgan Cambareri1, Alexis J. Yagielski1, Hanzhi T. Zhao1, Daniel A. Handwerker2,4, Javier Gonzalez-Castillo2,4, Peter A. Bandettini2,4, Elizabeth M. C. Hillman1,3,6,* Cell Reports (2023): https://doi.org/10.1016/j.celrep.2023.112527  
  
1. Zuckerman Mind Brain Behavior Institute, Department of Biomedical Engineering, Columbia University, New York, NY, USA  
2. Section on Functional Imaging Methods, Laboratory of Brain and Cognition, National Institute of Mental Health, Bethesda, MD, USA  
3. Department of Radiology, Columbia University Irving Medical Center, New York, NY, USA  
4. Functional MRI Core Facility, National Institutes of Health, Bethesda, MD, USA  
5. These authors contributed equally  
6. Lead contact  
*Correspondence: elizabeth.hillman@columbia.edu  
  
  PI: Elizabeth M. C. Hillman, PhD  
  Herbert and Florence Irving Professor  
  Mortimer B. Zuckerman Mind Brain Behavior Institute,  
  Professor of Biomedical Engineering and Radiology,  
  Columbia University,  
  Jerome L. Greene Science Center,  
  3227 Broadway, L5 Quad 5B,  
  New York, NY 10027  
  Lab Phone: 212-853-1097  
    
    The data used for this code are deposited at:http://dx.doi.org/10.17632/xd93nswg6h.1.
  
  
  **CODE OVERVIEW**  
  This GitHub repository consists of two main folders: Preprocessing and Analysis. Here is a brief overview of each folder and its contents:  
  
  _Preprocessing Folder_:  
  This folder contains a script that outlines the fundamental steps for preprocessing raw WFOM data. This code can be run on an example run by downloading 'cm128_day8_runB_all_data_shared.mat' and 'Spectra_share.mat' from http://dx.doi.org/10.17632/xd93nswg6h.1. This repository also includes behavioral camera recordings from the same run. 
  
  _Analysis Folder_:  
  Within this folder, you can find scripts related to conducting correlation analysis and non-negative least squares fits on the preprocessed WFOM data. The analysis is divided into five consecutive steps, where each step relies on the results of the previous step(s). It indicates that the analysis is iterative and builds upon the outcomes of earlier stages. To ensure the successful execution of these analysis steps, it is important to add the "Auxiliary_Code" folder to the MATLAB path. This code can be run on the WFOM extracted signals from http://dx.doi.org/10.17632/xd93nswg6h.1.
