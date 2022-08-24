#ifndef GENERAL_H_
#define GENERAL_H_

class InitializeInfo{
    public:
        InitializeInfo(DM dmcell, DM dmvertx);

    private:
        int         dim_;
        vector<int> localsize_;
        vector<int> localstart_;
        vector<int> globalsize_;
        vector<int> ghost_cell_;
        vector<int> ghost_vertx_;

}

#endif
