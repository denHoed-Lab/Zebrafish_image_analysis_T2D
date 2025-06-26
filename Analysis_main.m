%
%
% VAST image analysis pipeline for InsH2Bfabp background.
% 
% This pipeline was developed for the study:
% "In vivo, large-scale perturbation screen of T2D candidate genes highlights potential therapeutic targets"
%
% Authors: Endrina Mujica1*, Hanqing Zhang1*, Anastasia Emmanouilidou1*, Naomi Cook1, 
% Christoph Metzendorf1, Eugenia Mazzaferro1, Manoj Bandaru1, Joao Campos Costa1, 
% Ghazal Alavioon1, Tiffany Klingström1, Chrysoula Zalamitsou1, Natalie van Zuydam1, 
% Tessa Dronkers1, Mehran Hariri1, Bobby Mathews1, Adrianna Pakula1, Antje K Rottner2, 
% Han Sun3, Joshua W. Knowles4,5, Erik Ingelsson1, Mark I. McCarthy6, Martijn van de Bunt7, 
% Klaus Stensgaard Frederiksen8, Anna L Gloyn2,3,5,9, Hannah L. Brooke10, Anders Larsson11, 
% Sara Gry Vienberg8, Jason Flannick12, Amin Allalou13,14,15*, Marcel den Hoed1*#
%
% * denotes equal contributions
% # denotes corresponding author
%
% Affiliations:
% 1. The Beijer Laboratory and Department of Immunology, Genetics and Pathology; Uppsala University and SciLifeLab; Uppsala; Sweden
% 2. Oxford Centre for Diabetes, Endocrinology & Metabolism; University of Oxford; Oxford; UK
% 3. Division of Endocrinology; Department of Pediatrics; Stanford School of Medicine; Stanford University; Stanford; CA; USA
% 4. Department of Medicine; Division of Cardiovascular Medicine and Cardiovascular Institute; Stanford University; Stanford; CA, USA
% 5. Stanford Diabetes Research Center; Stanford, CA; USA
% 6. Wellcome Centre for Human Genetics; University of Oxford; Oxford, UK
% 7. Cytoki Pharma; Søborg; Denmark
% 8. Novo Nordisk A/S; Måløv; Denmark
% 9. Department of Genetics; Stanford School of Medicine; Stanford, CA; USA
% 10. Medical Epidemiology; Department of Surgical Sciences; Uppsala University; Uppsala; Sweden
% 11. Department of Medical Sciences; Division of Clinical Chemistry; Uppsala University; Uppsala; Sweden
% 12. Metabolism Program; The Broad Institute of Harvard and MIT; Cambridge; MA; USA
% 13. Department of Information Technology; Division of visual information and interaction; Uppsala University; Uppsala; Sweden
% 14. SciLifeLab BioImage Informatics Facility; Uppsala University; Uppsala, Sweden
% 15. DanioReadout; Department of Immunology, Genetics and Pathology; Uppsala University; Uppsala, Sweden
%
% Abstract: Genome-wide association studies (GWAS) have identified >1000 type 2 diabetes (T2D)-associated loci, 
% but functional characterization of candidate genes is lagging. We used a three-stage experimental perturbation 
% approach in zebrafish larvae to pinpoint high confidence causal genes amongst 43 human candidate genes. 
% In stage 1 we identified 21 genes affecting ≥1 of five T2D traits when zebrafish orthologues were perturbed 
% in vivo using CRISPR/Cas9 (beta cell mass, beta cell average/total nuclear ins expression, liver fat content, 
% random glucose content). Genes affecting ≥2 T2D traits in zebrafish larvae were enriched (p=4.9E-2) for rare 
% variant associations with T2D in humans. Nine of 21 genes had little or no prior evidence for a role in T2D. 
% In stage 2, mutagenesis of three of these genes (sirt1, atp2a3 and poldip2) also affected fasting glucose content 
% in zebrafish larvae. In stage 3 we show that effects of mutations in sirt1 and poldip2 are likely not driven by 
% global effects on early embryonic development. Mutations in sirt1 also tended to slow down clearance of glucose 
% after a 1 h glucose challenge. Together, CRISPR/Cas9 gene perturbation, semi-automated live imaging and deep 
% neural networks for image analysis yielded nine previously unanticipated T2D genes of which at least two are 
% likely to directly act on the pancreatic islet. Our approach improves our understanding of T2D etiology and 
% can yield new targets for intervention.
%
% Key words: CRISPR/Cas9, zebrafish larvae, fluorescence microscopy, deep learning
%
% Corresponding author: Marcel den Hoed (marcel.den_hoed@igp.uu.se)
%
% Version 1.1, 20250619
%
% Authors for this code: 
% Hanqing Zhang (hanzha@kth.se), 
% Amin Allalou (amin.allalou@it.uu.se), 
% Endrina Mujica (endrina.mujica@igp.uu.se), 
% Marcel den Hoed (marcel.den_hoed@igp.uu.se)
%
clear all; close all;

%% Inputs Pathes
source_path{1}='.\Test data\Test_experiment_demo';

%% load dependencies
if(~isdeployed)
    cPath=mfilename('fullpath');
    [GUIpath,~]=fileparts(cPath);
    lPath_f={'functions','_BetaCells','_Liver','_VAST','Config'};
    for lPath_i=lPath_f
        lPath_i=lPath_i{1};
        lPath_i=strcat(GUIpath,'\',lPath_i);
        addpath(genpath(lPath_i))
    end
end
%% Load Configuration
Config = Config_Demo();

%% Analysis
for fID=1:numel(source_path)
    inputFolder = source_path{fID};
    % Betacell Analysis
    BetalogFile=betaCells_framework(inputFolder,Config);
    % Liver Analysis
    LiverlogFile=Liver_framework(inputFolder,Config);
    % Vast image Analysis
    VAST_process(inputFolder,Config);
end
