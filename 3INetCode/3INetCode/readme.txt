# Deep learning assisted reconstruction of triangle-beam interference SIM images

The training, validation and testing data of four typical subcellular structures are available at the Zenodo link provided in the manuscript. 

To train or test the model, please first download the data, and place the data in the Data folder. Training widefield and GT data should be placed in Data\Raw and Data\GT, respectively; Validation widefield and GT data should be placed in Data\ValRaw and Data\ValGT, respectively; Testing widefield and GT data should be placed in Data\TestRaw and Data\TestGT, respectively. 

To train a new model, please run the Training.py file

To test the model, please run the Test.py file

If you want to use our pretrained model, please download it in the Zenodo link provided in the manuscript, and please the "final_network.pth" file in the Pretrained folder, and run the Test.py file, by alternating the default dir for model.load_state_dict function