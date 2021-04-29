#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <limits>
#include "use_hmm.h"
#include "model.h"

using namespace hmm;

void Model::rand_init_param() {

    double val = -log(n);
    // fill A: sum_j a_ij = 1

    srand(2);
    for (int i=0; i<n; i++)
    {
        double sum = 0.0;
        for (int j=0; j<n ; j++)
        {
            A[i][j] =  rand() % 10 + 1; // random number between 1 and 10
            sum += A[i][j];
        }
        for (int j=0; j<n ; j++)
        {
            A[i][j] =  log(A[i][j]) - log(sum);
        }
    }

    // fill B: \sum_k b[j][k] = 1
    val = -log(M);
    for (int j=0; j<n; j++)
        for (int k=0; k<M; k++)
            B[j][k]=val;

    // fill π: sum_i = 1
    val = -log(n);
    for (int i=0; i<n; i++)
        pi[i]=val;
}

void Model::init_param_EOS() 
{
    if (dic_obs->count(EOS) > 0)
    {
        // P(O_t=EOS|q_t=EOS) = 1
        for (int i=0; i<M; i++)
            B[(*dic_state)[EOS]][i] = -100.0;
        B[(*dic_state)[EOS]][(*dic_obs)[EOS]] = 1.0;
    }
}

double Model::ln_of_sum(double* ln_val, int start, int end) {
    // assertion: end>start
    double res;
    if (end-start > 1)
    {
        res = ln_val[start] + log((double)1.0 + (double)exp(ln_of_sum(ln_val, start+1, end) - ln_val[start]));
    } else {
        // end-start=1
        res = ln_val[start] + log((double)1.0 + (double)exp(ln_val[end] - ln_val[start]));
    }
    return res;
}

double Model::ln_of_sum(double ln_val1, double ln_val2) {
    return ln_val1 + log((double)1.0 + (double)exp(ln_val2 - ln_val1));
}

double** Model::alpha(vector<int>* seq_obs)
{  
    double test = 0.0;
    double *pt1 = &test;
    double **pt2 = &pt1;
    return pt2;
}

double** Model::beta(vector<int>* seq_obs){
    double test = 0.0;
    double *pt1 = &test;
    double **pt2 = &pt1;
    return pt2;
}

vector<int>* Model::viterbi(vector<int>* seq_obs)
{
    unsigned T = seq_obs->size();
    double delta[T][n];
    unsigned psi[T][n];

    cout << " shape : " << sizeof(delta)/sizeof(delta[0]) << " " << sizeof(delta[0])/sizeof(delta[0][0]) << endl;

    for(unsigned i = 0; i < n; i++)
    {   
        delta[0][i] = pi[i] + B[i][seq_obs->at(0)];
        psi[0][i] = 0;
    }

    for(unsigned t = 1; t < T; t++)
    {
        for(unsigned j = 0; j <= n - 1; j++)
        {
            vector<double>* tmp = new  vector<double>();
            for(unsigned x = 0; x <= n - 1; x++)
                tmp->push_back(delta[t - 1][x] + A[x][j]);

            delta[t][j] = get_max(tmp) + B[j][seq_obs->at(t)];
            psi[t][j] = argmax(tmp);
        }
    }

    vector<int>* q_star = new vector<int>();
    vector<double>* tmp = new vector<double>();

    for(unsigned i = 0; i < n; i++)
        tmp->push_back(delta[T-1][i]);
    
    double max = get_max(tmp);
    q_star->insert(q_star->begin(), argmax(tmp));

    for (int t = T - 2; t >= 0; t--)
        q_star->insert(q_star->begin(), psi[t + 1][q_star->at(0)]);

    /*cout << "Observation :" << endl;
    for(const int &val: *seq_obs)
         cout << val << " ";
    cout << endl;

    cout << "Météo :" << endl;
    for(const int &val: *q_star)
         cout << val << " ";
    cout << endl;*/

    return q_star;
}

double Model::get_max(vector<double>* list)
{
    double max = -9999999;
    for(const double &val : *list)
        if(val > max)
            max = val;

    return max;
}

unsigned Model::argmax(vector<double>* list)
{
    double max = -9999999;
    unsigned cpt = 0;
    for(const double &val : *list)
    {
        if(val > max)
        {
            max = val;            
            cpt++;
        }
    }
    cpt = cpt - 1;
    return cpt;
}

vector<double>* Model::posterior(vector<int>* seq_obs, vector<int>* seq_state) 
{
    vector<double>* test;
    return test;
}

vector<vector<int>*>* Model::predict(vector<vector<int>*>* data_obs) 
{
    vector<vector<int>*>* res = new vector<vector<int>*>();
    for (int k=0; k<data_obs->size(); k++)    
        res->push_back(Model::viterbi((*data_obs)[k]));
    return res;
}

double Model::pgen(vector<int>* seq_obs, double** alpha) 
{
    return 0.0;
}

double*** Model::xi(vector<int>* seq_obs, double* likelihood) 
{
    // assert that seq_obs is not empty
    double*** res = new double**[seq_obs->size()-1];
    for (int t=0; t<seq_obs->size()-1; t++)
    {
        res[t] = new double*[n];
        for (int i=0; i<n; i++)
            res[t][i] = new double[n];
    }

    double** alpha = Model::alpha(seq_obs);
    double** beta = Model::beta(seq_obs);

    // compute P(seq_obs)
    double denom = Model::pgen(seq_obs, alpha);

    for (int t=0; t<seq_obs->size()-1; t++)
    {
        for (int i=0; i<n; i++)
            for (int j=0; j<n; j++)
            {
                res[t][i][j] = alpha[t][i] + A[i][j] + B[j][(*seq_obs)[t+1]] + beta[t+1][j] - denom;
              }
        delete alpha[t];
        delete beta[t];
    }
    delete alpha[seq_obs->size()-1];
    delete alpha;
    delete beta[seq_obs->size()-1];
    delete beta;

    *likelihood = denom; 
    return res;
}

double Model::baum_welch_EMstep(vector<vector<int>*>* data)
{
    double log_likelihood_data = 0.0;

    double** new_A_num = new double*[n];
    double* new_AB_denom = new double[n];
    for(int i = 0; i < n; i++)
    {
        new_A_num[i] = new double[n];
        new_AB_denom[i] = 0.0;
        for (int j=0; j<n; j++)
        {
            new_A_num[i][j] = 0.0;
        }
    }
    double** new_B_num = new double*[n];
    for(int i = 0; i < n; i++)
    {
        new_B_num[i] = new double[M];
        for (int j=0; j<M; j++)
        {
            new_B_num[i][j] = 0.0;
        }
    }
    double* new_pi_num = new double[n];
    for (int i=0; i<n; i++)
        new_pi_num[i] = 0.0;

    for (int k=0; k<data->size(); k++)
    {
        vector<int>* seq_obs = (*data)[k];
        double likelihood;

        /* --------------- E step ------------------ */
        double*** xi = Model::xi(seq_obs, &likelihood);
        log_likelihood_data = ln_of_sum(log_likelihood_data,likelihood);
       
        // compute gamma
        double** gamma = new double*[seq_obs->size()-1];
        for (int t=0; t<seq_obs->size()-1; t++)
        {
            gamma[t] = new double[n];
            for (int i=0; i<n; i++)
                gamma[t][i] = ln_of_sum(xi[t][i],0,n-1);
        }

        /* --------------- M step ------------------ */
        for (int i=0; i<n; i++)
        {
            new_pi_num[i] =  ln_of_sum(new_pi_num[i], gamma[0][i]);
            double sum[seq_obs->size()-1];
            for (int t=0; t<seq_obs->size()-1; t++)
                    sum[t] = gamma[t][i];
            new_AB_denom[i] = ln_of_sum(new_AB_denom[i], ln_of_sum(sum,0,seq_obs->size()-2));
                
            for (int j=0; j<n; j++)
            {
                double sum2[seq_obs->size()-1];
                for (int t=0; t<seq_obs->size()-1; t++)
                    sum2[t] = xi[t][i][j];
                new_A_num[i][j] = ln_of_sum(new_A_num[i][j], ln_of_sum(sum2,0,seq_obs->size()-2));
            }
        }

        for (int j=0; j<n; j++)
        {
            for (int t=0; t<seq_obs->size()-1; t++)
                new_B_num[j][(*seq_obs)[t]] = ln_of_sum(new_B_num[j][(*seq_obs)[t]], gamma[t][j]);
        }


        for (int t=0; t<seq_obs->size()-1; t++)
        {
            for (int i=0; i<n; i++)
                delete xi[t][i];
            delete xi[t];
            delete gamma[t];
        }
        delete xi;
        delete gamma;
    }
    /* --------------- Update ------------------ */
    bool stop = true;
    for (int i=0; i<n; i++)
    {   
        double denom = new_AB_denom[i];
        for (int j=0; j<n ; j++)
        {
            double val = new_A_num[i][j] - denom;
            A[i][j] = val;
        }
        for (int k=0; k<M; k++)
        {
            double val = new_B_num[i][k] - denom;
            B[i][k] = val;
        }
        double val = new_pi_num[i] - (double)log(data->size());
        pi[i] = val;
    }
    
    /* --------------- Clean data -------------- */
    for(int i = 0; i < n; i++)
    {
        delete new_A_num[i];
        delete new_B_num[i];
    }
    delete new_A_num;
    delete new_B_num;
    delete new_AB_denom;
    delete new_pi_num;

    return log_likelihood_data;
}

void Model::train_baum_welch(vector<vector<int>*>* data, string filename) 
{
    rand_init_param();
    init_param_EOS();
    int num_iter = 0;
    save_hmm(filename);

    cout << "saving logs in baum-welch.log..." << endl;
    ofstream filout("baum-welch.log");
    bool stop = false;
    double prev_log_likelihood = Model::baum_welch_EMstep(data);
    if (filout)
            filout << num_iter << ";" << prev_log_likelihood << endl;
    while(!stop)
    {
        double cur_log_likelihood = Model::baum_welch_EMstep(data);
        num_iter++;        
        if (filout)
            filout << num_iter << ";" << cur_log_likelihood << endl;
        stop = abs(cur_log_likelihood - prev_log_likelihood) < STOP_CRITERION;
        prev_log_likelihood = cur_log_likelihood;

        if (num_iter%200==0)
        {
            save_hmm(filename+"_"+to_string(num_iter));
            cout << "(" << num_iter << ")" << flush;
        }
        else if (num_iter%100==0)
            cout << "*" << flush;
    }
    cout << endl;
    filout.close();
    init_param_EOS();
    save_hmm(filename);
}

bool Model::IsAlreadyExist(vector<int>* seq, int current)
{
    for(const int &val: *seq)
    {
        if(val == current)
            return true;
    }
    return false;
}

void Model::train_MLE(vector<vector<int>*>* data_obs, vector<vector<int>*>* data_state, string filename)
{
    /* --------------- Learn pi -------------- */
    // count starting states
    //int count_start_state[n];

    // calculer pi
    for(unsigned t = 0; t < n; t++)
    {
        double cpt_total = 0;
        double cpt_motif = 0;

        for(const vector<int>* vec : *data_state)
        {
            for(const int &val : *vec)
            {
                cpt_total++;
                if(val == t)
                    cpt_motif++;
            }
        }
        double result = cpt_motif/cpt_total;
        pi[t] = result;
        //cout << pi[t] << endl;
    }

    // calculer B
    //cout<< data_state->size() << " " << data_obs->size() << endl;
    for(unsigned r = 0; r < M; r++)
    {
        for(unsigned c = 0; c < n; c++)
        {
            double cpt_total = 0;
            double cpt_motif = 0;

            unsigned i = 0;
            for(const vector<int>* vec : *data_obs)
            {
                vector<int>* state = data_state->at(i);
                unsigned j = 0;
                for(const int &val : *vec)
                {
                    if(state->at(j) == c)
                        cpt_total++;

                    if(val == r && state->at(j) == c)
                        cpt_motif++;
                    j++;
                }
                i++;
            }
            //cout << "motif : " <<  cpt_motif << " total : " << cpt_total << endl;
            double result;
            if(cpt_motif == 0)
                result = 0;
            else
                result = cpt_motif/cpt_total;
            B[c][r] = result;
            //cout << B[c][r] << endl;
        }
    }
    
    // calculer A
    for(unsigned r = 0; r < n; r++)
    {
        for(unsigned c = 0; c < n; c++)
        {
            double cpt_total = 0;
            double cpt_motif = 0;

            unsigned i = 0;
            for(const vector<int>* vec : *data_state)
            {
                unsigned j = 0;
                for(const int &val : *vec)
                {
                    vector<int>* state = data_state->at(i);

                    if(j < state->size() - 1)
                    {
                        if(val == r && state->at(j + 1) == c)
                            cpt_motif++;
                        if(val == r)
                            cpt_total++;
                    }
                    j++;
                }
                i++;
            }
            //cout << "motif : " <<  cpt_motif << " total : " << cpt_total << endl;
            double result;
            if(cpt_motif == 0)
                result = 0;
            else
                result = cpt_motif/cpt_total;
            A[r][c] = result;
            //cout << A[r][c] << endl;
        }
    }
    /* --------------- Learn A -------------- */
    // count bigrams s_{i}s_{j} and unigrams s_{i}
    //int count_bigram_state[n][n];
    //int count_unigram_state[n];

    /* --------------- Learn B -------------- */
    // count pairs s_{i}o_{j}
    //int count_state_obs[n][M];

    /*for(unsigned i = 0;  i < n; i++)
        cout << pi[i] << " ";*/

    init_param_EOS();
    save_hmm(filename);
}


Model* Model::load_hmm(string filename) 
{
    Model* res;
    ifstream filin(filename);
    if (filin.is_open())
    {
        string line;
        int n = -1;
        int M = -1;
        if (getline (filin,line))
            n = stoi(line.substr(line.find("n=")+2));
        if (getline (filin,line))
            M = stoi(line.substr(line.find("M=")+2));

        if (n != -1 && M != -1)
        {
            res = new Model(n,M);
            
            enum Part {Part_start=0, Part_A=1, Part_B=2, Part_pi=3, Part_states=4, Part_obs=5};
            Part section = Part_start;
            int i = 0;
            while ( getline (filin,line) )
            {
                if (!line.empty()) {
                    if (line.rfind("# *** A",0) == 0)
                    {
                        section = Part_A;
                        i=0;
                    }    
                    else if (line.rfind("# *** B",0) == 0)
                    {
                        section = Part_B;
                        i=0;
                    }
                    else if (line.rfind("# *** pi",0) == 0)
                        section = Part_pi;
                    else if (line.rfind("# *** state",0) == 0)
                    {
                        section = Part_states;
                        i=0;
                    }
                    else if (line.rfind("# *** observation",0) == 0)
                    {
                        section = Part_obs;
                        i=0;
                    } 
                    else
                    {
                        if (section == Part_states) 
                        {
                            (*res->dic_state)[line] = i;
                            res->tab_state[i] = line;
                            i++;
                        }
                        else if (section == Part_obs)
                        {
                            (*res->dic_obs)[line] = i;
                            res->tab_obs[i] = line;
                            i++;
                        }
                        else 
                        {
                            int j = 0;
                            int maxJ = n;
                            if (section == Part_B)
                                maxJ = M;
                            istringstream iss(line);
                            string value;
                            while (getline(iss, value, ' ') && j<maxJ) 
                            {
                                // Replace , by . in case local settings saved ,
                                // for float numbers
                                int index = value.find(',');
                                if (index != string::npos)
                                    value.replace(index, 1, ".");
                                
                                if (section == Part_A || section == Part_B || section == Part_pi)
                                {
                                    double v = stod(value);
                                    double lv = log(v);
                                    if (v==0)
                                        lv = LOG_MINUS_INFTY;
                                    if (section == Part_A)
                                        res->A[i][j] = lv;
                                    else if (section == Part_B)
                                        res->B[i][j] = lv;
                                    else if (section == Part_pi)
                                        res->pi[j] = lv;
                                }
                                j++;
                            }
                            i++;
                        }
                    }
                }
            }
        }                
    }
    return res;
}

void Model::save_hmm(string filename) {
    ofstream filout(filename);
    if (filout.is_open())
    {
        filout << "n=" << n << endl;
        filout << "M=" << M << endl;
        filout << "# *** A parameters ***" << endl;
        for (int i=0; i < n; i++)
        {
            for (int j=0; j< n; j++)
                filout << A[i][j] << " ";
            filout << endl;
        }
        filout << "# *** B parameters ***" << endl;
        for (int i=0; i < n; i++)
        {
            for (int j=0; j< M; j++)
                filout << B[i][j] << " ";
            filout << endl;
        }
        filout << "# *** pi parameters ***" << endl;
        for (int i=0; i < n; i++)
        {
            filout << pi[i] << " ";
        }
        filout << endl;
        filout << "# *** state vocabulary ***" << endl;
        for (int i=0; i<n; i++)
            filout << tab_state[i] << endl;

        filout << "# *** observation vocabulary ***" << endl;
        for (int i=0; i<M; i++)
            filout << tab_obs[i] << endl;

        filout.close();
    }
}
