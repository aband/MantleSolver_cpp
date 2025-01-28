#ifndef TRANS_PARAM_H_
#define TRANS_PARAM_H_

int advfunc(const vector<double>& uin, const vector<double>& uout,
                  vector<double>& fin,       vector<double>& fout,
            const vector<vertex>& param,     vector<double>& LF,
            const vertex& unitnormal);

#endif
