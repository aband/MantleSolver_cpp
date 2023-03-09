            //! Uniformally add reconstruction methods.
            void AddReconstMethod(vector<indice> brm) {baseReconstMethod_.push_back(brm); AddWgts_();}; 
            void AddReconstMethod(vector<vector<indice>> brms) 
                                 {for (int i=0; i<brms.size(); i++){
                                      AddReconstMethod(brms[i]);     
                                  }}; 

            //! Two ways of updateing non linear weights with updated smoothness indicators.
            void UpdateNonLinearWgts(const MeshInfo& mi);
            void UpdateTwoStageNonLinearWgts(const MeshInfo& mi);

            // Collective methods
            void ResetWgts_();
            void AddWgts_();

            void UpdateNonLinearWgts_(const MeshInfo& mi, indice start);
            void UpdateFirstStageNonLinearWgts_(const MeshInfo& mi, indice start);
            void UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, indice start);

            vector< singleLevelReconstruction *> allLevels_;
            vector<vector<indice>> baseReconstMethod_; 
            vector< map<int, double> > linearWgts_;

            map<int, vector< map<int, double>>* > allLinearWgts_;
            map<int, vector< map<int, double> > > nonLinearWgts_;




void multiLevelReconstruction::AddWgts_() {

    map<int, double> lw;
    map<int, int> bias;

    const int sizeX = allLevels_[baseReconstMethod_.size()-1]->GetSizeX();

    for (auto& b: baseReconstMethod_[baseReconstMethod_.size()-1]){
             
        lw.insert({FlatIndic(sizeX,b),1.0});
        bias.insert({FlatIndic(sizeX,b),0});
    }

    linearWgts_.push_back(lw);
    etaBias_.push_back(bias);
}

void multiLevelReconstruction::ResetWgts_() {
    // Initialize linear weights and non linear weights
    // with the given information on reconstruction method
    linearWgts_.clear();
    etaBias_.clear();

    linearWgts_.resize(baseReconstMethod_.size());
    etaBias_.resize(baseReconstMethod_.size());

    for (int i=0; i<baseReconstMethod_.size();i++){
        const int sizeX = allLevels_[i]->GetSizeX();
        for (auto& b : baseReconstMethod_[i]){
            linearWgts_[i].insert({FlatIndic(sizeX, b),1.0});
            etaBias_[i].insert({FlatIndic(sizeX,b),0});
        }
    }

}

void multiLevelReconstruction::UpdateNonLinearWgts_(const MeshInfo& mi, indice start){

    assert(baseReconstMethod_.size() == allLevels_.size());

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                double value = linearWgts_[l].at(FlatIndic(sizeX,i))/ 
                               //pow(sm + scale*scale*eps0_ , max(sizeX, sizeY)) * 
                               pow(sm + scale*scale*eps0_ , sizeX+sizeY) * 
                               pow(eps0_*scale / sm+eps0_*
                               scale, etaBias_[l].at(FlatIndic(sizeX,i)));
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});

}

void multiLevelReconstruction::UpdateFirstStageNonLinearWgts_(const MeshInfo& mi, indice start){

    assert(baseReconstMethod_.size() == allLevels_.size());

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                int power = 0;
                if (allLevels_[l]->GetSizeX() * allLevels_[l]->GetSizeY() == 1){
                    power = 1;} else {
                    power = 2;
                }

                double value = linearWgts_[l].at(FlatIndic(sizeX,i))/ 
                               pow(sm + scale*scale*eps0_ , power); 
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});

}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts_(const MeshInfo& mi, indice start){

    UpdateFirstStageNonLinearWgts_(mi,start);

    vector< map<int,double> > nlw(linearWgts_.size());

    double sum = 0.0;

    for (int l=0; l<allLevels_.size(); l++){
       const int sizeX = allLevels_[l]->GetSizeX();
       const int sizeY = allLevels_[l]->GetSizeY();
       for (auto& i:baseReconstMethod_[l]){
            indice owner = start + i;
            if (allLevels_[l]->CheckExist(mi, owner)){
                double scale = allLevels_[l]->GetScale(FlatIndic(mi,owner));
                double sm = allLevels_[l]->CalculateSmoothnessIndic(mi,owner);
                // Get updated smoothness indicators
                vector< map<int, double> > oldnlw = nonLinearWgts_.at(FlatIndic(mi,start));
                double value = oldnlw[l].at(FlatIndic(sizeX,i))/ 
                               //pow(sm + scale*scale*eps0_ , max(sizeX, sizeY)) * 
                               pow(sm + scale*scale*eps0_ , sizeX+sizeY) * 
                               pow(eps0_*scale / sm+eps0_*
                               scale, etaBias_[l].at(FlatIndic(sizeX,i)));
                nlw[l].insert({FlatIndic(sizeX,i) , value});
                sum += value;
            }            

        }
    }

    for (int l=0; l<allLevels_.size(); l++){
        if (nlw[l].empty() ==0){
            for (auto & in:nlw[l]){
                in.second = in.second/sum; 
            }
        }
    }

    nonLinearWgts_.erase(FlatIndic(mi,start));
    nonLinearWgts_.insert({FlatIndic(mi,start) , nlw});
}

void multiLevelReconstruction::UpdateNonLinearWgts(const MeshInfo& mi){
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice add {i,j};
        indice start = mi.MPIlocalCellStart+add;
        UpdateNonLinearWgts_(mi,start);
    } }
}

void multiLevelReconstruction::UpdateTwoStageNonLinearWgts(const MeshInfo& mi){
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
    for (int i=0; i<mi.MPIlocalCellSize[0]; i++){
        indice add {i,j};
        indice start = mi.MPIlocalCellStart+add;
        UpdateTwoStageNonLinearWgts_(mi,start);
    } }
}


