#include <htool/htool.hpp>


using namespace std;
using namespace htool;


class Mycondition final : public VirtualAdmissibilityCondition {
  public:
    bool ComputeAdmissibility(const VirtualCluster &target, const VirtualCluster &source, double eta) const override {
        bool admissible = 2 * std::min(target.get_rad(), source.get_rad()) < eta * std::max((norm2(target.get_ctr() - source.get_ctr()) - target.get_rad() - source.get_rad()), 0.);
        return admissible;
    }
};
void produitP(Block<double>* B, const vector<double>& x, vector<double>& y, char op){
    int of_t = B->get_target_cluster().get_offset();
    int of_s = B->get_source_cluster().get_offset();
    int sz_t = B->get_target_cluster().get_size();
    int sz_s = B->get_source_cluster().get_size();
    if (!(B->get_block_data() ==nullptr) and (B->nb_sons()==0)){
        if (op =='N'){
            B->get_block_data()->add_mvprod_row_major(x.data()+of_s,y.data()+of_t,1,'N',op);
        }
        else{ B->get_block_data()->add_mvprod_row_major(x.data()+of_t,y.data()+of_s,1,'N',op);
        }
    }
    else{
        for (int r =0; r<B->nb_sons();++r){
            Block<double>& Br = B->get_son(r);
            produitP(&Br,x,y,op);
        } }}

 class SumExpression {
  private :
    vector<Matrix<double>> SR;
    vector<int> off;
    vector<Block<double>*> SH;

  public:
    //Constructeur qui prend deux blocs
    SumExpression(Block<double>* H, Block<double>* K){ vector<Block<double>*> vb; vb.push_back(H); vb .push_back(K); SH = vb;}
    // Constructeur avec tout
    SumExpression(const vector<Matrix<double>> SR0,const vector<int> off0,const vector<Block<double>*> Sh0){this->SR=SR0;this->off=off0;this->SH=Sh0;}

    // fonction restrict
    SumExpression Restrict(const VirtualCluster& t , const VirtualCluster& s){
        //cout<<"on entre dans restrict"<< endl;
        int oft = t.get_offset(); int ofs = s.get_offset(); int szt = t.get_size(); int szs = s.get_size();
        //cout<< oft<<','<<ofs<<','<<szt<<','<<szs<<endl;
        vector<Matrix<double>> SR0 = this->SR;
        vector<int> off0 = this-> off;
        vector<Block<double>*> SH0 = this-> SH;
        vector<Matrix<double>> sr;vector<int> of;
        // taille des deux (on a deux matrices a chaques fois
        int nr = SR0.size()/2; int nh = SH0.size()/2;
        // réstriction de SR;
        for(int k = 0; k < nr; ++k){
            Matrix<double> U = SR0[2*k]; Matrix<double> V = SR0[2*k+1];
            int ofu = off0[2*k]; int ofv = off0[2*k+1];
            if( (oft >=ofu) and ( szt <= U.nb_rows()) and (ofs >= ofv) and ( szt <= V.nb_cols())){
                    // si on  a des blocs concernés on extrait
                    // ------> implémenter un extract rapide
                    Matrix<double> Uk; Matrix<double> Vk;
                    for (int i = 0 ; i< szt; ++i){
                        for(int j =0; j < U.nb_cols(); ++j){
                            Uk(i,j) = U(oft-ofu +i, ofs-ofv+j);
                        }
                    }
                    for (int i = 0 ; i< V.nb_cols(); ++i){
                        for(int j =0; j < szs; ++j){
                            Vk(i,j) = V(oft-ofu +i, ofs-ofv+j);
                        }
                    }
                    sr.push_back(Uk);sr.push_back(Vk); of.push_back(oft); of.push_back(ofs);}
                }

            //La on a la restriction de Sr
            // //on veut Sh maintentnant
            // On parcours les HK de SH
            vector<Block<double>*> vh;
            vector<Block<double>*> vk;
            vector<Block<double>*> vtemp;
            for (int k =0; k< nh; ++k){


                Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
                //On prend rho=r = H.source= K.target();
                int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
                int nb_r = K->get_target_cluster().get_nb_sons();
                int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
                int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
                //On regarde si ils sont pas déja a la bonne taille

                if (((oft== of_h) and ( szt== sz_h)) and ( (ofs==of_k) and (szs==sz_k)))
                {
                    vtemp.push_back(H); vtemp.push_back(K);}
                else{
                    // on regarde les concernés
                    //cout<< H->nb_sons()<<","<<K->nb_sons()<<endl;
                    for (int i = 0; i< H->nb_sons();++i){
                        Block<double>& Hk = H->get_son(i);
                        for(int j=0; j< K->nb_sons(); ++j){
                            Block<double>& Kk = K->get_son(j);
                            if ( ((oft>=Hk.get_target_cluster().get_offset() ) and (szt<= Hk.get_target_cluster().get_size())) and ((ofs>=Kk.get_source_cluster().get_offset() ) and (szs<= Kk.get_source_cluster().get_size()))){
                                if((Hk.get_source_cluster().get_offset()==Kk.get_target_cluster().get_offset()) and (Hk.get_source_cluster().get_size()==Kk.get_target_cluster().get_size())){
                                    vtemp.push_back(&Hk); vtemp.push_back(&Kk); }	} }

                    }} }
            //normalement vtemp c'est notre nouveaux Sh mais il a peut être des low rank
            //On parcours et on regarde si il en a a des low rank;
            int nn = vtemp.size()/2;
            //cout<<"!!!!!!!!!!!!!!!!!!"<<nn<<endl;
            for (int l = 0; l< nn ; ++l){
                Block<double>* H = vtemp[2*l]; Block<double>* K= vtemp[2*l+1];
            // on regarde si les deux sont low rank
                if(!(H->get_low_rank_block_data() ==nullptr) and !(K->get_low_rank_block_data() == nullptr)){
                    Matrix<double> Uh,Vh,Uk,Vk;
                    Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
                    Uk = K->get_low_rank_block_data()->Get_U();Uk = K->get_low_rank_block_data()->Get_V();
                    Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
                    int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
                    sr.push_back(U); sr.push_back(V); of.push_back(rt); of.push_back(rs);
                }
                //Celle de droite low rank
                else if( !(K->get_low_rank_block_data() == nullptr)){
                    Matrix<double> U = K->get_low_rank_block_data()->Get_U();
                    Matrix<double> V = K->get_low_rank_block_data()->Get_V();
                    Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
                    for (int rr = 0 ; rr <U.nb_cols();++rr){
                        const vector<double> x = U.get_col(rr);
                        int rp = H->get_target_cluster().get_size();
                        vector<double> y (rr,0.);
                        produitP(H,x,y,'N');
                        for(int kl = 0 ; kl< y.size(); ++kl){
                            W(kl,rr)=y[kl];}}
                    sr.push_back(W);sr.push_back(V); of.push_back(H->get_target_cluster().get_offset()); of.push_back(K->get_source_cluster().get_offset());}
                //celle de guauche low rank
                else if( !(H->get_low_rank_block_data()==nullptr)){
                    Matrix<double> U = H->get_low_rank_block_data()->Get_U();
                    Matrix<double> V = H->get_low_rank_block_data()->Get_V();
                    Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());
                    for (int rr = 0 ; rr <V.nb_rows();++rr){
                        const vector<double> x = V.get_row(rr);
                        vector<double> y (K->get_target_cluster().get_size(),0.);
                        produitP(K,x,y,'T');
                        for(int kl = 0 ; kl< y.size(); ++kl){
                            W(rr,kl)=y[kl];}
                    }
                    sr.push_back(U);sr.push_back(W); of.push_back(H->get_target_cluster().get_offset()); of.push_back(K->get_source_cluster().get_offset());
                }
                else{
                    //ni H ni K ne sont low rank on les laisse dans vh
                    vh.push_back(H);vh.push_back(K);}
            }


            SumExpression Res(sr,of,vh);
            //cout<<"on sort de retrict"<< endl;
            return Res;}
//getters
    vector<Matrix<double>> get_sr(){return this->SR;}
    vector<int> get_off(){return this->off;}
    vector<Block<double>*> get_sh(){return this->SH;}

    //fonction evaluate
    // on l'appelle qua quand la feuille est pas adm donc a priori on est pas low rank
    // trois cas : cd serait bloc*bloc(bizzard car ce erait pas le meme critére de compression
    //ou bien cloc*plein =plein *bloc = iteration(bloc vecteur)

    Matrix<double> eval_sr(){
        vector<Matrix<double>> sr = this->SR;
        Matrix<double> res;
        for(int k =0; k< sr.size()/2; ++k){
            Matrix<double> U= sr[2*k];Matrix<double> V = sr[2*k+1];
            res= res + U*V;}
        return res;
    }

    Matrix<double> eval_sh(){
        vector<Block<double>*> sh  =this->SH;
        Matrix<double> val( sh[0]->get_target_cluster().get_size(),sh[1]->get_source_cluster().get_size());
        for (int k =0; k < sh.size()/2; ++k){
            Block<double>* H =sh[2*k];Block<double>* K = sh[2*k+1];
            if ((H->nb_sons()==0) or (K->nb_sons()==0)){
                if (!(K->get_dense_block_data() == nullptr)) {
                    cout << "ici" << endl;
                    const Matrix<double> *ktemp = K->get_dense_block_data();
                    Matrix<double> kkk          = *ktemp;
                    // H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                    for (int l = 0; l < kkk.nb_cols(); ++l) {
                        const vector<double> xk = kkk.get_col(l);
                        vector<double> yk(H->get_target_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'N');
                        for (int ii = 0; ii < H->get_target_cluster().get_size(); ++ii) {
                            val(ii , l ) += yk[ii];
                        }
                    }
                }
                // si celle de droite a un virtual block
                // du coup on fait un produit ligne block
                else {
                    cout << "la" << endl;
                    const Matrix<double> *ktemp = H->get_dense_block_data();
                    Matrix<double> kkk = *ktemp;
                    for (int l = 0; l < kkk.nb_rows(); ++l) {
                        const vector<double> xk = kkk.get_row(l);
                        vector<double> yk(K->get_source_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'T');
                        cout << "prod ok" << endl;
                        for (int ii = 0; ii < K->get_source_cluster().get_size(); ++ii) {
                            val(l, ii) += yk[ii];
                        }
                    }
            }

        }
        else{
            for(int nt = 0 ; nt< H->get_target_cluster().get_nb_sons(); ++nt){
                for(int ns = 0; ns< K->get_source_cluster().get_nb_sons(); ++ns){
                    Matrix<double> valk (H->get_target_cluster().get_son(nt).get_size(),K->get_source_cluster().get_son(ns).get_size());
                    SumExpression sk = this->Restrict(H->get_target_cluster().get_son(nt),K->get_source_cluster().get_son(ns));
                    valk = sk.eval_sh();
                    for(int i = 0 ; i< valk.nb_rows();++i){
                        for(int j =0; j< valk.nb_cols();++j){
                            val(i+nt*(H->get_target_cluster().get_son(1).get_offset()-H->get_target_cluster().get_son(0).get_offset()),j+ns*(K->get_source_cluster().get_son(1).get_offset()-H->get_source_cluster().get_son(0).get_offset()))
                                += valk(i,j);

                        }
                    }
                }
            }
                 }
        }
        return val;
        }

    Matrix<double> eval(){return this->eval_sh()+this->eval_sr();}

    Matrix<double> evaluate(int oft,int ofs){
        vector<Matrix<double>> sr = this->SR;
        vector<int> of = this->off;
        vector<Block<double>*> sh = this->SH;
        Matrix<double> val;
        for( int k =0; k< sr.size()/2;++k){
            val= val + sr[2*k]*sr[2*k+1];
        }
        for(int k = 0 ; k < sh.size(); ++k){
        Block<double>* H = sh[2*k];Block<double>* K = sh[2*k+1];
        // si au moins un est pas subd
        if( (H->nb_sons()==0 )or (K->nb_sons()==0)){
                if (!(K->get_dense_block_data() == nullptr)) {
                    cout << "ici" << endl;
                    const Matrix<double> *ktemp = K->get_dense_block_data();
                    Matrix<double> kkk          = *ktemp;
                    // H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                    for (int l = 0; l < kkk.nb_cols(); ++l) {
                        const vector<double> xk = kkk.get_col(l);
                        vector<double> yk(H->get_target_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'N');
                        for (int ii = 0; ii < H->get_target_cluster().get_size(); ++ii) {
                            val(ii + H->get_target_cluster().get_offset()-oft, l + K->get_source_cluster().get_offset()-ofs) += yk[ii];
                        }
                    }
                }
                // si celle de droite a un virtual block
                // du coup on fait un produit ligne block
                else {
                    cout << "la" << endl;
                    const Matrix<double> *ktemp = H->get_dense_block_data();
                    Matrix<double> kkk = *ktemp;
                    for (int l = 0; l < kkk.nb_rows(); ++l) {
                        const vector<double> xk = kkk.get_row(l);
                        vector<double> yk(K->get_source_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'T');
                        cout << "prod ok" << endl;
                        for (int ii = 0; ii < K->get_source_cluster().get_size(); ++ii) {
                            val(l+H->get_target_cluster().get_offset()-oft, ii+K->get_source_cluster().get_offset()-ofs) += yk[ii];
                        }
                    }
                    cout << "assemblage ok" << endl;
                    //val = val + temp;
                }
            }
        else{
            for (int i = 0 ; i < H->nb_sons(); ++i){
                for(int j = 0; j< K->nb_sons(); ++j){

                }
            }
        }
        }
        return val;
    }

    void Evaluate(Matrix<double> M,int oft, int ofs) {
        int nr                     = M.nb_rows();
        int nc                     = M.nb_cols();
        vector<Matrix<double>> sr  = this->SR;
        vector<int> of             = this->off;
        vector<Block<double> *> sh = this->SH;
        Matrix<double> val;
        // on parcours les lowrank
        for (int k = 0; k < sr.size() / 2; ++k) {
            Matrix<double> U = sr[2 * k];
            Matrix<double> V = sr[2 * k + 1];
            val              = val + U * V;
        }

        // on parcours sh
        double *vh;
        for (int k = 0; k < sh.size() / 2; ++k) {
            Block<double> *H = sh[2 * k];
            Block<double> *K = sh[2 * k + 1];

            Matrix<double> temp(H->get_target_cluster().get_size(), K->get_source_cluster().get_size());
            // on vérifie qu'ils soit pas tout les deux des blocs
            cout<<"!!!!!"<<H->nb_sons()<<','<<K->nb_sons()<< endl;
            if ((H->nb_sons() == 0 )or (K->nb_sons() == 0)) {
                // normalement avec mvmprod si on met mu=2 on a pas besoin de différencier les cas
                // fut quand meme savoir le quel prend le mvprod
                // si les les deux sont pleines ont fait blas
                // Si celle de gauche a un virtual block
                if (!(K->get_dense_block_data() == nullptr)) {
                    cout << "ici" << endl;
                    const Matrix<double> *ktemp = K->get_dense_block_data();
                    Matrix<double> kkk          = *ktemp;
                    // H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                    for (int l = 0; l < kkk.nb_cols(); ++l) {
                        const vector<double> xk = kkk.get_col(l);
                        vector<double> yk(H->get_target_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'N');
                        for (int ii = 0; ii < H->get_target_cluster().get_size(); ++ii) {
                            M(ii + H->get_target_cluster().get_offset()-oft, l + K->get_source_cluster().get_offset()-ofs) += yk[ii];
                        }
                    }
                }
                // si celle de droite a un virtual block
                // du coup on fait un produit ligne block
                else {
                    cout << "la" << endl;
                    const Matrix<double> *ktemp = H->get_dense_block_data();
                    Matrix<double> kkk = *ktemp;
                    for (int l = 0; l < kkk.nb_rows(); ++l) {
                        const vector<double> xk = kkk.get_row(l);
                        vector<double> yk(K->get_source_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'T');
                        cout << "prod ok" << endl;
                        for (int ii = 0; ii < K->get_source_cluster().get_size(); ++ii) {
                            M(l+H->get_target_cluster().get_offset()-oft, ii+K->get_source_cluster().get_offset()-ofs) += yk[ii];
                        }
                    }
                    cout << "assemblage ok" << endl;
                    //val = val + temp;
                }
            }
            else{
                cout<<"recursion"<<endl;
                SumExpression& S =*this;
                for(int nt = 0; nt < H->get_target_cluster().get_nb_sons(); ++nt){
                    for (int ns = 0 ; ns < K->get_source_cluster().get_nb_sons(); ++ns){
                        cout<<"k,l"<<nt<<','<<ns<<endl;
                        cout<<"on restreint a "<< endl;
                        cout<< H->get_target_cluster().get_nb_sons()<<','<<K->get_source_cluster().get_nb_sons()<<endl;
                        SumExpression Sk = S.Restrict(H->get_target_cluster().get_son(nt),K->get_source_cluster().get_son(ns));
                        cout<<"restrict ok"<< endl;
                        Sk.Evaluate(M,oft,ofs);
                    }
                }
            }
        }
    }

/*
    Matrix<double> Evaluate(int oft, int ofs){
        vector<Matrix<double>> sr = this->SR;
        vector<int> of = this->off;
        vector<Block<double>*> sh = this->SH;
        Matrix<double> val;
        //on parcours les lowrank
        for(int k =0; k < sr.size()/2;++k){
            Matrix<double> U= sr[2*k]; Matrix<double> V =sr[2*k+1];
            val= val +U*V;

        }
        cout<< "low rank ok"<<endl;
        cout<< sh.size()/2<< endl;

        //on parcours sh
        double* vh;
        for(int k = 0; k < sh.size()/2;++k){
            cout<<"k:"<<k<<endl;
            Block<double>* H = sh[2*k];Block<double>* K = sh[2*k+1];
            cout<< H->get_target_cluster().get_size()<<','<<K->get_source_cluster().get_size()<< endl;

            Matrix<double> temp ( H->get_target_cluster().get_size(),K->get_source_cluster().get_size());
            // normalement avec mvmprod si on met mu=2 on a pas besoin de différencier les cas
            //fut quand meme savoir le quel prend le mvprod
            //si les les deux sont pleines ont fait blas
            // Si celle de gauche a un virtual block
            if( (H->nb_sons()==0) or (K->nb_sons()==0)){
            if ( !(K->get_dense_block_data()==nullptr)){
                cout<< "ici"<<endl;
                //K->get_block_data()->add_mvprod_row_major(H->get_dense_block_data(),temp,2,'N','N');
                const Matrix<double>* ktemp = K->get_dense_block_data();
                Matrix<double> kkk = *ktemp;
                //H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                for (int l = 0; l< kkk.nb_cols();++l){
                    const vector<double> xk = kkk.get_col(l);
                    vector<double> yk(H->get_target_cluster().get_size(),0);
                    produitP(H,xk,yk,'N');
                    for (int ii = 0 ; ii < H->get_target_cluster().get_size();++ii) {
                        temp(ii, l) = yk[ii];
                    }
                }
                val = val+temp;
            }
            // si celle de droite a un virtual block
            // du coup on fait un produit ligne block
            else{
                cout<<"la"<<endl;
                const Matrix<double>* ktemp = H->get_dense_block_data();
                cout<< H->nb_sons()<<','<<K->nb_sons()<<endl;
                cout<< ktemp->nb_rows()<<endl;
                Matrix<double> kkk = *ktemp;
                cout<<"?"<< endl;
                //H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                for (int l = 0; l< kkk.nb_rows();++l){
                    cout<<"l:"<<l<< endl;
                    const vector<double> xk = kkk.get_row(l);
                    vector<double> yk(K->get_source_cluster().get_size(),0);
                    produitP(H,xk,yk,'T');
                    cout<<"prod ok"<< endl;
                    for (int ii = 0 ; ii < K->get_source_cluster().get_size();++ii){
                        temp(l,ii) = yk[ii];
                    }
                }
                cout<< "assemblage ok"<<endl;
                val = val+temp;}
            //K->get_block_data()->add_mvprod_row_major(H->get_dense_block_data(),temp,2,'N','T');
            cout<<"k:"<<k<<endl;
        }
        else{
            for (int st = 0; st< H->get_target_cluster().get_nb_sons();++st){
                for(int ss = 0 ; ss < K->get_source_cluster().get_nb_sons();++ss){
                    SumExpression Sk  =  this->Restrict(H->get_target_cluster().get_son(st),K->get_source_cluster().get_son(ss));
                    Sk.Evaluate(oft,ofs);
                }
            }
        }
        }
        //Matrix<double> temp ( H->get_target_cluster().get_size(),K->get_source_cluster().get_size());
        //H.add_mvprod_row_major(K->get_)


        return val;
    }
*/
};
class MyMatrix : public VirtualGenerator<double> {
    const vector<double> &p1;
    const vector<double> &p2;
    int space_dim;

  public:
    // Constructor
    MyMatrix(int space_dim0, int nr, int nc, const vector<double> &p10, const vector<double> &p20) : VirtualGenerator(nr, nc), p1(p10), p2(p20), space_dim(space_dim0) {}

    // Virtual function to overload
    double get_coef(const int &k, const int &j) const {
        return (1.) / (4 * M_PI * std::sqrt(1e-5 + std::inner_product(p1.begin() + space_dim * k, p1.begin() + space_dim * k + space_dim, p2.begin() + space_dim * j, double(0), std::plus<double>(), [](double u, double v) { return (u - v) * (u - v); })));
    }

    // Virtual function to overload
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const override {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                ptr[j + M * k] = this->get_coef(rows[j], cols[k]);
            }
        }
    }

    // Matrix vector product
    std::vector<double> operator*(std::vector<double> a) {
        std::vector<double> result(nr, 0);
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                result[j] += this->get_coef(j, k) * a[k];
            }
        }
        return result;
    }

    // Frobenius norm
    double norm() {
        double norm = 0;
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                norm += this->get_coef(j, k);
            }
        }
        return norm;
    }


/*
  //////////////////////////////////////////////////
 //  std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> restrict(std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> S, Cluster t, Cluster s){
//   std::vector<Block<T>> Sr = S.first;
//   //il faut faire la restriction des low rank a t,s
//   std::vector<Block<T>> Sr0;    int offset_t = t.get_offset(); int offset_s = s.get_offset(); int size_t = t.get_size(); int size_s = s.get_size();
//   for (auto it= Sr.begin(); it!= Sr.end(); ++it){
//     Block<T> lr = *it;
//     int repere = 0;
//     // ____________________//
//     //Pas sure que cette méthode marche mais sinon il faut mettre la matrice en full, faire sa restriction et refaire ACA;
//     // on peut aussi juste faire une boucle for k < lr.get_nb_son()
//     //____________________//
//     // On récupère la réstriction a tau sigma en vérifiant sur chaques blocs des get_son() si leurs clusters sont les mêmes
//     while(!(&lr.get_son(repere) ==nullptr)){
//       Block<T> &lr_son = lr.get_son(repere);
//       int offset_t0=lr_son.get_target_cluster().get_offset(); int offset_s0= lr_son.get_source_cluster.get_offset();
//       int size_t0= lr_son.get_target_cluster().get_size(); int size_s0 = lr_son.get_source_cluster().get_size();
//       // on vérifie que les cluster sont les mêmes <=> même sze et offset
//       if( ((offset_t == offset_t0) and (offset_s == offset_s0))  and ((size_t0 == size_t) and (size_s0 == size_s))){
// 	Sr0.push_back(lr_son);}
//       repere +=1;
//     }
//   }
//   //on initialise Sh0 a 0
//   std::vector<pair<Block<T>,Block<T>>> Sh = S.second;
//   std::vector<pair<Block<T>,Block<T>>> Sh0;
//   //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
//   for (auto it= Sh.begin(); it!= Sh.end(); ++it){
//     pair<Block<T>,Block<T>> HK = *it; Block<T> H=HK.first; Block<T> K = HK.second;
//     Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
//     int repere1 = 0;
//     //on parcours les enfants de rho=r
//     while(!(r.get_son_ptr(repere1) == nullptr)){
//       //on cherche les restrictions de H et K parmi leur fils comme pour Sr
//       int repereH =0;int repereK=0;Block<T> H0;Block<T> K0;
//       int offset_r0 = r.get_son_ptr(repere).get_offset();
//       int size_r0= r.get_son_ptr(repere).get_size();
//       //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
//       while(!(&H.get_son(repereH) == nullptr)){
// 	Block<T> &Htemp = H.get_son(repereH);
// 	int offset_t0 = Htemp.get_target_cluster().get_offset();
// 	int offset_s0 = Htemp.get_source_cluster().get_offset();
// 	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
// 	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
// 	  H0 = Htemp;}
// 	repereH +=1;}
//       // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
//       while(!(&K.get_son(repereK) == nullptr)){
// 	Block<T> &Ktemp = K.get_son(repereK);
// 	int offset_t0 = Ktemp.get_target_cluster().get_offset();
// 	int offset_s0 = Ktemp.get_source_cluster().get_offset();
// 	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
// 	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
// 	  K0 = Ktemp;}
// 	repereK +=1;}
//       //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
//       // Si c'est le cas le produit est forcément low rank
// 	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
// 	//On veut faire ABt = H0*K0------>push_back->Sr0
// 	// trois possibilité: lr*lr, lr*hmat, hmat*lr
// 	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
// 	// Dans les trois cas on a un lr en résultat
// 	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
// 	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
// 	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
// 	 LowRankMatrix<T>& hr0 = H0.get_low_rank_block_data();
// 	 LowRankMatrix<T>& kr0 = K0.get_low_rank_block_data();
// 	 DenseBlockData<T>& hf0 = H0.get_dense_block_data();
// 	 DenseBlockData<T>& kf0 = K0.get_dense_block_data();
// 	//cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
// 	if((hf0 == nullptr) and (kf0 == nullptr)){
// 	  LowRankMatrix<T> hk0 = h0*k0;}
// 	//cas 2 : hmat*lowrank
// 	else if( hr0 == nullptr ){
// 	  Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
// 	  int rank = U.nb_cols();
// 	  LowRankMatrix<T> hk0;
// 	  Matrix<T>* hk0U = hk0.Get_U();
// 	  Matrix<T>* hk0V = res.Get_V();
// 	  *hk0V = V;
// 	  Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
// 	  for ( int k =0; k< rank;++k){
// 	    vector<T> Uk = U.get_col(k);
// 	    vector<T> temp = H0*Uk;
// 	    for ( int l =0; l< U.nb_rows();++l){
// 	      mattemp(l,k)= Uk[l];
// 	    }
// 	  }
// 	  *hk0U = mattemp;
// 	}
// 	//cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
// 	// Sinon on code un produit vecteur mat xA = At x
// 	else if( kr0 ==nullptr){
	  
// 	}
	  
//       }
//       else{
// 	pair<Block<T>,Block<T>> HK0; HK0.first = H0;HK0.second = K0;
// 	Sh0.push_back(HK0);
//       }
//     }
//   }
//   std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> S0;
//   S0.first = Sr0; S0.second = Sh0;
//   return S0;
// }  
*/
  ////////////////////////////////////



  
};



//produit entre deux matrices pleines


Matrix<double> Mprod(Matrix<double> A, Matrix<double> B){
    Matrix<double> res(A.nb_rows(),B.nb_cols());
    for (int l = 0; l < B.nb_cols();++l){
        vector<double> x = B.get_col(l); vector<double> y = A*x;
        for(int k =0; k< A.nb_rows(); ++k){
            res(k,l) = y[k];
        }
    }
    return res;}
// Produit de prierre --------------------> évidemment ca marche mais attention il faut faire:
//global_to_cluster(&Rt->get_source_cluster(),xx.data(),xx_perm.data());
//cluster_to_global(&Rt->get_target_cluster(),yy0.data(),yy_perm.data());
//produit(Rt,xx_perm,yy_perm);
//cluster_to_global(&Rt->get_target_cluster(),yy_perm.data(),yy0.data());

void produit(Block<double>* B, const vector<double>& x, vector<double>& y){
    int of_t = B->get_target_cluster().get_offset();
    int of_s = B->get_source_cluster().get_offset();
    int sz_t = B->get_target_cluster().get_size();
    int sz_s = B->get_source_cluster().get_size();
    if (!(B->get_block_data() ==nullptr) and (B->nb_sons()==0)){
        B->get_block_data()->add_mvprod_row_major(x.data()+of_s,y.data()+of_t,1,'N','N');
    }
    else{
        for (int r =0; r<B->nb_sons();++r){
            Block<double>& Br = B->get_son(r);
            produit(&Br,x,y);
        } }}

/*
// Produit transpo avec la méthode de prierre
void produitP(Block<double>* B, const vector<double>& x, vector<double>& y, char op){
    int of_t = B->get_target_cluster().get_offset();
    int of_s = B->get_source_cluster().get_offset();
    int sz_t = B->get_target_cluster().get_size();
    int sz_s = B->get_source_cluster().get_size();
    if (!(B->get_block_data() ==nullptr) and (B->nb_sons()==0)){
        if (op =='N'){
        B->get_block_data()->add_mvprod_row_major(x.data()+of_s,y.data()+of_t,1,'N',op);
    }
    else{ B->get_block_data()->add_mvprod_row_major(x.data()+of_t,y.data()+of_s,1,'N',op);
    }
    }
    else{
        for (int r =0; r<B->nb_sons();++r){
            Block<double>& Br = B->get_son(r);
            produitP(&Br,x,y,op);
        } }}

 */
 // Produit transpo
// Produit marche mais j'ia quand même une petite erreur ...
void Produit(Block<double>* B, const vector<double> x, vector<double>* y,int* n){
  int of_t = B->get_target_cluster().get_offset(); int of_s = B->get_source_cluster().get_offset();
  int sz_t = B->get_target_cluster().get_size(); int sz_s = B->get_source_cluster().get_size();
  // if( (B->get_block_data()==nullptr) and(B->nb_sons()==0)){cout<<"!!!!!!!!!"<<endl;}
  //if ((B->nb_sons()==0) and (B->get_block_data()== nullptr)){cout << *n<< endl;}
  //if (B->nb_sons()==0){
  if (!(B->get_block_data() ==nullptr) and (B->nb_sons()==0)){
      int r = *n; r= r+1; *n=r;
    //if( !(B->get_block_data()==nullptr)){
    vector<double> xt;
    for (int k =0; k< sz_s;++k){
      xt.push_back(x[k+of_s]);}
    vector<double> yy(sz_t,0);
    vector<double> py= *y;
    // for (int kk =0; kk< yy.size();++kk){
    //   yy[kk]=py[kk];}
    B->get_block_data()->add_mvprod_row_major(x.data(),yy.data(),1,'N','N');
    for (int k =0 ; k< sz_t; ++k){
      py[k+of_t]+=yy[k];}
    *y = py;
  }
  //}
  else{
   for (int r =0; r<B->nb_sons();++r){
     Block<double>& Br = B->get_son(r);
     vector<double> yy (y->size(),0);
     vector<double>* yz = &yy;
     Produit(&Br,x,y,n);
     // vector<double> py= *y;
     // vector<double> temp = *yz;
     // for (int k =0 ; k< Br.get_target_cluster().get_size(); ++k){
     //   py[k+Br.get_target_cluster().get_offset()]+=temp[k];}
     // *y=py;
     } }}
// Produit transpo
void Produitt(Block<double>* B, const vector<double> x, vector<double>* y,int* n){
  int of_t = B->get_target_cluster().get_offset(); int of_s = B->get_source_cluster().get_offset();
  int sz_t = B->get_target_cluster().get_size(); int sz_s = B->get_source_cluster().get_size();
  // if( (B->get_block_data()==nullptr) and(B->nb_sons()==0)){cout<<"!!!!!!!!!"<<endl;}
  //if ((B->nb_sons()==0) and (B->get_block_data()== nullptr)){cout << *n<< endl;}
  //if (B->nb_sons()==0){
  if (!(B->get_block_data() ==nullptr) and (B->nb_sons()==0)){
      int r = *n; r= r+1; *n=r;
    //if( !(B->get_block_data()==nullptr)){
    vector<double> xt;
    for (int k =0; k< sz_t;++k){
      xt.push_back(x[k+of_t]);}
    vector<double> py= *y;
    for (int k0=0; k0 < sz_s; ++ k0){
      double tempp =0;
      // C'est la qu'il faut changer
      for (int k1 =0; k1 <sz_t; ++k1){
	vector<double> yy(sz_t,0);
	vector<double> xtk(sz_s,0);xtk[k0]=xt[k1]; 
	B->get_block_data()->add_mvprod_row_major(xtk.data(),yy.data(),1,'N','N');
	tempp += yy[k1];}
      py[k0+of_t]+= tempp;}
    *y = py;
  }
  //}
  else{
   for (int r =0; r<B->nb_sons();++r){
     Block<double>& Br = B->get_son(r);
     vector<double> yy (y->size(),0);
     vector<double>* yz = &yy;
     Produit(&Br,x,y,n);
     // vector<double> py= *y;
     // vector<double> temp = *yz;
     // for (int k =0 ; k< Br.get_target_cluster().get_size(); ++k){
     //   py[k+Br.get_target_cluster().get_offset()]+=temp[k];}
     // *y=py;
     } }}

// fonction restrict pour les low rank Ca normalement ca marche
//cette fonction est appelé uniquement avec des matrices
pair<vector<Matrix<double>>,vector<int>> restrict_lr(pair<vector<Matrix<double>>,vector<int>> SR, const VirtualCluster &t, const VirtualCluster &s){
  pair<vector<Matrix<double>>,vector<int>> Res;
  vector<Matrix<double>> sk; vector<int> ofk;
  
    int oft = t.get_offset(); int szt = t.get_size(); int ofs = s.get_offset(); int szs = s.get_size();
    vector<Matrix<double>> sr = SR.first; vector<int> of = SR.second;
    int n = sr.size()/2;
    //pair<pair<int, Matrix<double>>,pair<int , Matrix<double>>> repere = SR[k];
    //pair< int,Matrix<double>> Uref = repere.first; int ofu = Uref.first; Matrix<double> U = Uref.second; int szU = U.nb_rows();
    //pair< int,Matrix<double>> Vref = repere.second; int ofv = Vref.first; Matrix<double> V = Vref.second; int szV = V.nb_cols();
    for (int k = 0; k < n ; ++k){
      Matrix<double> U = sr[2*k]; Matrix<double> V =  sr[2*k+1];
      int ofu = of[2*k]; int ofv = of[2*k+1];
      if( ((oft >= ofu) and (szt <= U.nb_rows())) and  ((ofs >= ofv) and (szs <= V.nb_cols()))){
	// repere c'est pour savoir a partir d'ou on rargarde dan la grande matrice -> nouvel offset
	// offset de target pour U et source pour V
	  int repere_i = oft-ofu; int repere_j = ofs-ofv;
	  Matrix<double> Uk( szt,U.nb_cols()); Matrix<double> Vk(V.nb_rows(),szs);
	  for (int i = 0; i< szt; ++i){
	    for(int j =0; j < U.nb_cols(); ++j){
	      Uk(i,j)=U(i+repere_i,j);}
	  }
	  for(int i= 0; i < V.nb_rows(); ++i){
	    for(int j =0 ; j< szs; ++j){
	      Vk(i,j)= V(i,j+repere_j);}
	  }
	  //maintenant on a juste a push nos matrices réduites et nos nouveaux offset
	  sk.push_back(Uk);
	  sk.push_back(Vk);
	  ofk.push_back(ofu+repere_i);ofk.push_back(ofv+repere_j);}
  }
    SR.first = sr ; SR.second = of;

  return SR;}
/*
//fonction restrict
pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict(pair < vector<Matrix<double>>,vector<int>>SR0, vector<Block<double>*>SH0,const VirtualCluster& t, const VirtualCluster& s){
  int of_t = t.get_offset(); int sz_t = t.get_size(); int of_s = s.get_offset(); int sz_s = s.get_size();
  pair<vector<Matrix<double>>,vector<int>> Sr = restrict_lr(SR0,t,s);
  pair<vector<Matrix<double>>,vector<int>> SR;
  pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Res;
  Res.first = SR0;Res.second= SH0;
  vector<Matrix<double>> Suv = SR0.first; vector<int> off = SR0.second;
  int n = SH0.size()/2;
  int oft = t.get_offset(); int ofs = s.get_offset();
  int szt = t.get_size(); int szs = s.get_size();
  // On parcours les HK de SH
      vector<Block<double>*> vh;
    vector<Block<double>*> vk;
    vector<Block<double>*> vtemp;

  for (int k =0; k< n; ++k){
    Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
    //On prend rho=r = H.source= K.target();
    int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
    int nb_r = K->get_target_cluster().get_nb_sons();
    int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
    int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
    //déja faut extraire les paires concerné
    if (((oft>= of_h) and ( szt<= sz_h)) and ( (ofs>=of_k) and (szs<=sz_k))){

	// on regarde le concerné
	for (int i = 0; i< H->nb_sons();++i){
	  Block<double>& Hk = H->get_son(i);
	  for(int j=0; j< K->nb_sons(); ++j){
	    Block<double>& Kk = K->get_son(j);
	    if ( ((oft>=Hk.get_target_cluster().get_offset() ) and (szt<= Hk.get_target_cluster().get_size())) and ((ofs>=Kk.get_source_cluster().get_offset() ) and (szs<= Kk.get_source_cluster().get_size()))){
	      if ((Kk.get_target_cluster().get_offset()==Hk.get_source_cluster().get_offset() ) and (Kk.get_target_cluster().get_size()== Hk.get_source_cluster().get_size())){
	      vtemp.push_back(&Hk); vtemp.push_back(&Kk);} }	} }
    }}
      //normalement vtemp c'est notre nouveaux Sh mais il a peut être des low rank
      int nn = vtemp.size()/2;
      for (int l = 0; l< nn ; ++l){
	Block<double>* H = vtemp[2*l]; Block<double>* K= vtemp[2*l+1];
	// on regarde si les deux sont low rank
	if(!(H->get_low_rank_block_data() ==nullptr) and !(K->get_low_rank_block_data() == nullptr)){
	  Matrix<double> Uh,Vh,Uk,Vk;
	  Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
	  Uk = K->get_low_rank_block_data()->Get_U();Uk = K->get_low_rank_block_data()->Get_V();
	  Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
	  int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
	  Suv.push_back(U); Suv.push_back(V); off.push_back(rt); off.push_back(rs);
	}
	//Celle de droite low rank
	else if( !(K->get_low_rank_block_data() == nullptr)){
	  Matrix<double> U = K->get_low_rank_block_data()->Get_U();
	  Matrix<double> V = K->get_low_rank_block_data()->Get_V();
	  Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
	  for (int rr = 0 ; rr <U.nb_cols();++rr){
	    vector<double> x = U.get_col(rr);
	    int rp = H->get_target_cluster().get_size();
	    vector<double> y (rr,0.);
	    vector<double>* py = &y;
	    Produit(H,x,py,0);
	    vector<double> xux = *py;
	    for(int kl = 0 ; kl< xux.size(); ++kl){
	      W(kl,rr)=xux[kl];}}
	  Suv.push_back(W);Suv.push_back(V); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());}
	//celle de guauche low rank
	else if( !(H->get_low_rank_block_data()==nullptr)){
	    Matrix<double> U = H->get_low_rank_block_data()->Get_U();
	  Matrix<double> V = H->get_low_rank_block_data()->Get_V();
	  Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());
	  for (int rr = 0 ; rr <V.nb_rows();++rr){
	    vector<double> x = V.get_row(rr);
	    vector<double> y (K->get_target_cluster().get_size(),0.);
	    vector<double>* py = &y;
	    Produitt(K,x,py,0);
	    vector<double> xux = *py;
	    for(int kl = 0 ; kl< xux.size(); ++kl){
	      W(kl,rr)=xux[kl];}
	  }
	  Suv.push_back(U);Suv.push_back(W); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());
	}
	  else{
	    vh.push_back(H);vh.push_back(K);}}

	  SR0.first = Suv; SR0.second = off; Res.first =SR0; Res.second  = vh;


 return Res;}
*/


//fonction restrict
pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict(pair<pair < vector<Matrix<double>>,vector<int>>, vector<Block<double>*>>S,const VirtualCluster& t, const VirtualCluster& s){
    int of_t = t.get_offset(); int sz_t = t.get_size(); int of_s = s.get_offset(); int sz_s = s.get_size();
    pair<vector<Matrix<double>>,vector<int>> Sr = S.first;
    vector<Block<double>*> SH0 = S.second;
    pair<vector<Matrix<double>>,vector<int>> SR0 = restrict_lr(Sr,t,s);
    pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Res;
    Res.first = SR0;Res.second= SH0;
    vector<Matrix<double>> Suv = SR0.first; vector<int> off = SR0.second;
    int n = (SH0.size())/2;
    int oft = t.get_offset(); int ofs = s.get_offset();
    int szt = t.get_size(); int szs = s.get_size();
    // On parcours les HK de SH
    vector<Block<double>*> vh;
    vector<Block<double>*> vk;
    vector<Block<double>*> vtemp;
    cout<<"n"<<n;
    for (int k =0; k< n; ++k){
        Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
        //On prend rho=r = H.source= K.target();
        int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
        int nb_r = K->get_target_cluster().get_nb_sons();
        int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
        int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
        //déja faut extraire les paires concerné
        if (((oft>= of_h) and ( szt<= sz_h)) and ( (ofs>=of_k) and (szs<=sz_k))){

            // on regarde les concernés
            for (int i = 0; i< H->nb_sons();++i){
                Block<double>& Hk = H->get_son(i);
                for(int j=0; j< K->nb_sons(); ++j){
                    Block<double>& Kk = K->get_son(j);
                    if ( ((oft>=Hk.get_target_cluster().get_offset() ) and (szt<= Hk.get_target_cluster().get_size())) and ((ofs>=Kk.get_source_cluster().get_offset() ) and (szs<= Kk.get_source_cluster().get_size()))){
                        if ((Kk.get_target_cluster().get_offset()==Hk.get_source_cluster().get_offset() ) and (Kk.get_target_cluster().get_size()== Hk.get_source_cluster().get_size())){
                            vtemp.push_back(&Hk); vtemp.push_back(&Kk);} }	} }
        }
    }
    //normalement vtemp c'est notre nouveaux Sh mais il a peut être des low rank
    int nn = vtemp.size()/2;
    for (int l = 0; l< nn ; ++l){
        Block<double>* H = vtemp[2*l]; Block<double>* K= vtemp[2*l+1];
        // on regarde si les deux sont low rank
        if(!(H->get_low_rank_block_data() ==nullptr) and !(K->get_low_rank_block_data() == nullptr)){
            Matrix<double> Uh,Vh,Uk,Vk;
            Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
            Uk = K->get_low_rank_block_data()->Get_U();Uk = K->get_low_rank_block_data()->Get_V();
            Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
            int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
            Suv.push_back(U); Suv.push_back(V); off.push_back(rt); off.push_back(rs);
        }
        //Celle de droite low rank
        else if( !(K->get_low_rank_block_data() == nullptr)){
            Matrix<double> U = K->get_low_rank_block_data()->Get_U();
            Matrix<double> V = K->get_low_rank_block_data()->Get_V();
            Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
            for (int rr = 0 ; rr <U.nb_cols();++rr){
                vector<double> x = U.get_col(rr);
                int rp = H->get_target_cluster().get_size();
                vector<double> y (rr,0.);
                vector<double>* py = &y;
                Produit(H,x,py,0);
                vector<double> xux = *py;
                for(int kl = 0 ; kl< xux.size(); ++kl){
                    W(kl,rr)=xux[kl];}}
            Suv.push_back(W);Suv.push_back(V); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());}
        //celle de guauche low rank
        else if( !(H->get_low_rank_block_data()==nullptr)){
            Matrix<double> U = H->get_low_rank_block_data()->Get_U();
            Matrix<double> V = H->get_low_rank_block_data()->Get_V();
            Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());
            for (int rr = 0 ; rr <V.nb_rows();++rr){
                vector<double> x = V.get_row(rr);
                vector<double> y (K->get_target_cluster().get_size(),0.);
                vector<double>* py = &y;
                Produitt(K,x,py,0);
                vector<double> xux = *py;
                for(int kl = 0 ; kl< xux.size(); ++kl){
                    W(kl,rr)=xux[kl];}
            }
            Suv.push_back(U);Suv.push_back(W); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());
        }
        else{
            //ni H ni K ne sont low rank on les laisse dans vh
            vh.push_back(H);vh.push_back(K);}}

    pair<vector<Matrix<double>>,vector<int>> sr; sr.first = Suv; sr.second = off;
    Res.first = sr; Res.second = vh;

    return Res;}



//VRAI fonction restrict ---------------------------> Ca a l'aire de marcher
pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict0(pair<pair < vector<Matrix<double>>,vector<int>>, vector<Block<double>*>>S,const VirtualCluster& t, const VirtualCluster& s){
    int of_t = t.get_offset(); int sz_t = t.get_size(); int of_s = s.get_offset(); int sz_s = s.get_size();
    pair<vector<Matrix<double>>,vector<int>> Sr = S.first;
    vector<Block<double>*> SH0 = S.second;
    //on fait direct la restriction a de Sr
    pair<vector<Matrix<double>>,vector<int>> SR0 = restrict_lr(Sr,t,s);
    pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Res;
    Res.first = SR0;Res.second= SH0;
    vector<Matrix<double>> Suv = SR0.first; vector<int> off = SR0.second;
    int n = (SH0.size())/2;
    int oft = t.get_offset(); int ofs = s.get_offset();
    int szt = t.get_size(); int szs = s.get_size();
    // On parcours les HK de SH
    vector<Block<double>*> vh;
    vector<Block<double>*> vk;
    vector<Block<double>*> vtemp;
    for (int k =0; k< n; ++k){
        Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
        //On prend rho=r = H.source= K.target();
        int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
        int nb_r = K->get_target_cluster().get_nb_sons();
        int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
        int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
        //On regarde si ils sont pas déja a la bonne taille

        if (((oft== of_h) and ( szt== sz_h)) and ( (ofs==of_k) and (szs==sz_k)))
        {
            vtemp.push_back(H); vtemp.push_back(K);}
        else{
            // on regarde les concernés
            cout<< H->nb_sons()<<","<<K->nb_sons()<<endl;
            for (int i = 0; i< H->nb_sons();++i){
                Block<double>& Hk = H->get_son(i);
                    for(int j=0; j< K->nb_sons(); ++j){
                        Block<double>& Kk = K->get_son(j);
                        if ( ((oft>=Hk.get_target_cluster().get_offset() ) and (szt<= Hk.get_target_cluster().get_size())) and ((ofs>=Kk.get_source_cluster().get_offset() ) and (szs<= Kk.get_source_cluster().get_size()))){
                            if((Hk.get_source_cluster().get_offset()==Kk.get_target_cluster().get_offset()) and (Hk.get_source_cluster().get_size()==Kk.get_target_cluster().get_size())){
                            vtemp.push_back(&Hk); vtemp.push_back(&Kk); }	} }

    }} }
    //normalement vtemp c'est notre nouveaux Sh mais il a peut être des low rank
    int nn = vtemp.size()/2;
    for (int l = 0; l< nn ; ++l){
        Block<double>* H = vtemp[2*l]; Block<double>* K= vtemp[2*l+1];
        // on regarde si les deux sont low rank
        if(!(H->get_low_rank_block_data() ==nullptr) and !(K->get_low_rank_block_data() == nullptr)){
            Matrix<double> Uh,Vh,Uk,Vk;
            Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
            Uk = K->get_low_rank_block_data()->Get_U();Uk = K->get_low_rank_block_data()->Get_V();
            Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
            int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
            Suv.push_back(U); Suv.push_back(V); off.push_back(rt); off.push_back(rs);
        }
        //Celle de droite low rank
        else if( !(K->get_low_rank_block_data() == nullptr)){
            Matrix<double> U = K->get_low_rank_block_data()->Get_U();
            Matrix<double> V = K->get_low_rank_block_data()->Get_V();
            Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
            for (int rr = 0 ; rr <U.nb_cols();++rr){
                vector<double> x = U.get_col(rr);
                int rp = H->get_target_cluster().get_size();
                vector<double> y (rr,0.);
                vector<double>* py = &y;
                Produit(H,x,py,0);
                vector<double> xux = *py;
                for(int kl = 0 ; kl< xux.size(); ++kl){
                    W(kl,rr)=xux[kl];}}
            Suv.push_back(W);Suv.push_back(V); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());}
        //celle de guauche low rank
        else if( !(H->get_low_rank_block_data()==nullptr)){
            Matrix<double> U = H->get_low_rank_block_data()->Get_U();
            Matrix<double> V = H->get_low_rank_block_data()->Get_V();
            Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());
            for (int rr = 0 ; rr <V.nb_rows();++rr){
                vector<double> x = V.get_row(rr);
                vector<double> y (K->get_target_cluster().get_size(),0.);
                vector<double>* py = &y;
                Produitt(K,x,py,0);
                vector<double> xux = *py;
                for(int kl = 0 ; kl< xux.size(); ++kl){
                    W(kl,rr)=xux[kl];}
            }
            Suv.push_back(U);Suv.push_back(W); off.push_back(H->get_target_cluster().get_offset()); off.push_back(K->get_source_cluster().get_offset());
        }
        else{
            //ni H ni K ne sont low rank on les laisse dans vh
            vh.push_back(H);vh.push_back(K);}}

    pair<vector<Matrix<double>>,vector<int>> sr; sr.first = Suv; sr.second = off;
    Res.first = sr; Res.second = vh;

    return Res;}

//fonction de merde pour envoyer le coeff d'un virtual block
double coef(int i , int j ,int n , int m, const VirtualBlockData<double>* val){
    vector<double> x(m,0); vector<double> y (n,0);x[i]=1;
    val->add_mvprod_row_major(x.data(),y.data(),1,'N','N');
    return y[j];
}

//Fonction de merde pour avoir une collonne ( on va assembler nos matrices en faisant des produits matices vecteurs
vector<double> Col(int i  ,int n , int m, const VirtualBlockData<double>* val){
    vector<double> y;
    for (int j =0 ; j<n ; ++j){
        y.push_back(coef(i,j,n,m,val));
    }
    return y;
}

//fonction de merde pour avoir la matrice
Matrix<double> get_mat(int n , int m , const VirtualBlockData<double>* A,const VirtualBlockData<double>* B){
    Matrix<double> Res (n,m);
    for(int l =0; l<m ; ++l){
        vector<double> y = Col(l,n,m,B);
        for(int k =0; k< n; ++n){
            Res(k,l) = y[k];}
    }
    return Res;}
// Fonction multiplication
// en gros on regarde garde le bloc Lk et on le stocke en void dans L si on le calcule
// Ca marche pas en initialisant avec une H mat donc peut être si on regarde jute sa computed block
void Prod(HMatrix<double>* L,Block<double>* Lr,pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> S,const VirtualCluster& t,const  VirtualCluster& s){
  if (!(Lr->nb_sons() ==0)){
    for (int k =0; k< Lr->nb_sons(); ++k){
      Block<double>& Lk = Lr->get_son(k);
      pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> S0 = Restrict0(S,Lk.get_target_cluster(),Lk.get_source_cluster());
      Prod(L,&Lk,S0,Lk.get_target_cluster(),Lk.get_source_cluster());
    }
  }
  else{
    pair<vector<Matrix<double>>,vector<int>>SR = S.first; vector<Block<double>*> Sh= S.second;
    vector<Matrix<double>> Sr = SR.first; vector<int> off = SR.second;
    //Cas ou on calcule parce que c'est full (normalement du coup il y a rien dans Sr
    if (!(Lr->IsAdmissible()) ) {
      for(int k =0; k < Sh.size()/2; ++k){
	Block<double>* H = Sh[2*k]; Block<double>* K = Sh[2*k+1];
	if ( !(H->get_block_data()==nullptr) and !( K->get_block_data()==nullptr) ) {
	  // H*K ? -> est_ce que je dois juste définir une nouvelle mvprod??
	}
      }
    }
    //sinon tout est dans sr
     else{
        
      vector<vector<double>> uu,vv;
      for (int i =0 ; i< Sr.size()/2 ; ++i){
	vector<double> uuu;
	Matrix<double> U = Sr[2*i]; Matrix<double> V = Sr[2*i+1];
	for(int l = 0; l< U.nb_cols(); ++l){
	  vector<double> uk  = U.get_col(l); uu.push_back(uk);}
	for (int l =0; l< V.nb_rows(); ++l){
	  vector<double> vk = V.get_row(l); vv.push_back(vk);}
      }
      // on fabrique U et V;
      Matrix<double> U (uu[0].size(),uu.size());
      for (int k =0;k< uu[0].size();++k){
	for(int l =0; l<uu.size(); ++l){
	  U(k,l) = uu[l][k];} }
      Matrix<double> V(vv.size(),vv[0].size());
      for (int k =0; k< vv.size();++k){
	for (int l=0; l< vv[0].size(); ++l){
	  V(k,l) =  vv[k][l];}}
      Lr->get_low_rank_block_data()->Get_U() = U;
      Lr->get_low_rank_block_data()->Get_V() = V;}}
}

	
void testProd(Block<double>* L,Block<double>* Lr,pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> S,const VirtualCluster& t,const  VirtualCluster& s) {
    if (!(Lr->nb_sons() == 0)) {
        for (int k = 0; k < Lr->nb_sons(); ++k) {
            Block<double> &Lk                                                           = Lr->get_son(k);
            pair<pair<vector<Matrix<double>>, vector<int>>, vector<Block<double> *>> S0 = Restrict0(S, Lk.get_target_cluster(), Lk.get_source_cluster());
            testProd(L, &Lk, S0, Lk.get_target_cluster(), Lk.get_source_cluster());
        }
    } else {
        pair<vector<Matrix<double>>, vector<int>> SR = S.first;
        vector<Block<double> *> Sh                   = S.second;
        vector<Matrix<double>> Sr                    = SR.first;
        vector<int> off                              = SR.second;
        // Cas ou on calcule parce que c'est full (normalement du coup il y a rien dans Sr
        if (!(Lr->IsAdmissible())) {
            for (int k = 0; k < Sh.size() / 2; ++k) {
                Block<double> *H = Sh[2 * k];
                Block<double> *K = Sh[2 * k + 1];
                if (!(H->get_block_data() == nullptr) and !(K->get_block_data() == nullptr)) {
                    // H*K ? -> est_ce que je dois juste définir une nouvelle mvprod??
                    //est ce que je vais devoir faire des

                }
            }
        }
        // sinon tout est dans sr
        else {

            vector<vector<double>> uu, vv;
            for (int i = 0; i < Sr.size() / 2; ++i) {
                vector<double> uuu;
                Matrix<double> U = Sr[2 * i];
                Matrix<double> V = Sr[2 * i + 1];
                for (int l = 0; l < U.nb_cols(); ++l) {
                    vector<double> uk = U.get_col(l);
                    uu.push_back(uk);
                }
                for (int l = 0; l < V.nb_rows(); ++l) {
                    vector<double> vk = V.get_row(l);
                    vv.push_back(vk);
                }
            }
            // on fabrique U et V;
            Matrix<double> U(uu[0].size(), uu.size());
            for (int k = 0; k < uu[0].size(); ++k) {
                for (int l = 0; l < uu.size(); ++l) {
                    U(k, l) = uu[l][k];
                }
            }
            Matrix<double> V(vv.size(), vv[0].size());
            for (int k = 0; k < vv.size(); ++k) {
                for (int l = 0; l < vv[0].size(); ++l) {
                    V(k, l) = vv[k][l];
                }
            }
            Lr->get_low_rank_block_data()->Get_U() = U;
            Lr->get_low_rank_block_data()->Get_V() = V;
        }
    }
}

//fonction restrict2
// pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict2(pair < vector<Matrix<double>>,vector<int>>  SR0, vector<Block<double>*>SH0,const VirtualCluster& t, const VirtualCluster& s){
//   int of_t = t.get_offset(); int sz_t = t.get_size(); int of_s = s.get_offset(); int sz_s = s.get_size();
//   //on fait directement la restrictions des low rank
//   pair<vector<Matrix<double>>,vector<int>> Sr = restrict_lr(SR0,t,s);
//   pair<vector<Matrix<double>>,vector<int>> SR;
//   pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Res;
//   Res.first = SR0;Res.second= SH0;
//   // on initialise Res 
//   int n = SH0.size()/2;
//   // On parcours les HK de SH
//   for (int k =0; k< n; ++k){
//     Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
//     //On prend rho=r = H.source= K.target();
//     int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
//     int nb_r = K->get_target_cluster().get_nb_sons();
//     int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
//     int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
//     // On regarde si l'un des deux  low rank
//     if ((H->get_low_rank_block_data() == nullptr)) and (K->get_low_rank_block_data() == nullptr)){
//     //aucun low rank 
	
	
//     //Cas ou son.target== t
//    //  if ( ((sz_h == sz_t ) and (sz_k== sz_s)) and (( of_h  == of_t) and (of_k == of_s))){
//   //     for (int i = 0; i< nb_r; ++i){
//   // 	// On test si l'un des deux est une feuille.
//   // 	if( (H->IsAdmissible()) or ( K->IsAdmissible())){
	  
//   // 	  //trois cas
//   // 	  if( !(H->get_low_rank_block_data() == nullptr) and !( K->get_low_rank_block_data()==nullptr)){
//   // 	    // Les deux sont low rank
//   // 	    Matrix<double> Uh = H->get_low_rank_block_data()->Get_U(); Matrix<double> Vh = H->get_low_rank_block_data()->Get_V();
//   // 	    Matrix<double> Uk = K->get_low_rank_block_data()->Get_U(); Matrix<double> Vk = K->get_low_rank_block_data()->Get_V();
//   // 	    Matrix<double> u = Uh; Matrix<double> v= Vh*(Uk*Vk);
//   // 	    cout<<"!!!"<<endl;
//   // 	    cout<< Sr.first.size()<<endl;
//   // 	    Sr.first.push_back(u);Sr.first.push_back(v);Sr.second.push_back(of_t);Sr.second.push_back(of_s);
//   // 	    cout<< Sr.first.size()<<endl;

//   // 	   	}
//   // 	  else if( !(K->get_low_rank_block_data()==nullptr)){
//   // 	    Matrix<double> Uk = K->get_low_rank_block_data()->Get_U();Matrix<double> Vk = K->get_low_rank_block_data()->Get_V();
//   // 	    int nc = Uk.nb_cols();
//   // 	    Matrix<double> (sz_t,nc);
//   // 	    for (int r1 = 0; r1 < nc ; ++r1){
//   // 	      vector<double> x = Uk.get_col(r1);
//   // 	      //La il faut faire la multiplication hmat vecteur mais on a que un bloc !
//   // 	      //vector<double> y = 
//   // 	    }}
//   // 	}
//   //     }
//   //   }
//   //   else{
//   //     cout<<"chelou"<<endl;}

//   // }
//   //for (int k =0; k < n ; endl){
//   //Block<double>* H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
    
//   //}
//   return Res;}

int Hmult (Block<double>* L,Block<double>* Lk, SumExpression S) {


    if (!((Lk->get_target_cluster().get_nb_sons() == 0) and (Lk->get_source_cluster().get_nb_sons() == 0))) {
        int rep = 0;
        for (int k = 0; k < Lk->get_source_cluster().get_nb_sons(); ++k) {
            for (int l = 0; l < Lk->get_target_cluster().get_nb_sons(); ++l) {
                SumExpression sk = S.Restrict(Lk->get_target_cluster().get_son(l), Lk->get_source_cluster().get_son(k));
                Lk->build_son(Lk->get_target_cluster().get_son(l), Lk->get_source_cluster().get_son(k));
                cout<<" pas de pb build"<< endl;
                cout << Lk->nb_sons() << endl;
                Hmult(L, &Lk->get_son(rep), sk);
                cout<<"pas de pb a appeler hmult"<<endl;
                rep += 1;
            }
        }
    } else {
        return
/*
        cout << "evaluate" << endl;
        Matrix<double> val(Lk->get_target_cluster().get_size(), Lk->get_source_cluster().get_size());
        S.Evaluate(val, Lk->get_target_cluster().get_offset(), Lk->get_source_cluster().get_offset());
        // Matrix<double> val = S.Evaluate(Lk->get_target_cluster().get_offset(), Lk->get_source_cluster().get_offset());
        cout << "evalua ok" << endl;
        // Matrix<double>* vv = &val;
        const Matrix<double> *m = Lk->get_dense_block_data();
        m                       = &val;
        cout << "affectation ok" << endl;

        //  = &val;
*/
        /*
        cout<< "on évalue"<< endl;
        vector<Matrix<double>> sr0 = S.get_sr(); vector<Block<double>*> sh0 = S.get_sh();
        Matrix<double> Mtemp;
        if (sr0.size()==0) {
            int nr = sh0[0]->get_target_cluster().get_size();
            int nc = sh0[1]->get_source_cluster().get_size();
            Matrix<double> Mmtemp(nr, nc); Mtemp= Mmtemp;}
        else{int nr = sr0[0].nb_rows();
            int nc = sr0[1].nb_cols();
            Matrix<double> Mmtemp(nr, nc); Mtemp =Mmtemp;}

        for(int k =0; k < sh0.size()/2; ++k) {
            Block<double> *H = sh0[2 * k];
            Block<double> *K = sh0[2 * k + 1];
            if ((H->nb_sons() > 0) and (K->nb_sons() > 0)) {
                Hmult(Lk, Lk, SumExpression(H, K));
            } else {
                if (!(K->get_dense_block_data() == nullptr)) {
                    cout << "ici" << endl;
                    const Matrix<double> *ktemp = K->get_dense_block_data();
                    Matrix<double> kkk          = *ktemp;
                    // H->get_block_data()->add_mvprod_row_major(ktemp,temp,2,'N','N');
                    for (int l = 0; l < kkk.nb_cols(); ++l) {
                        const vector<double> xk = kkk.get_col(l);
                        vector<double> yk(H->get_target_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'N');
                        for (int ii = 0; ii < H->get_target_cluster().get_size(); ++ii) {
                            Mtemp(ii, l) += yk[ii];
                        }
                    }
                }
                // si celle de droite a un virtual block
                // du coup on fait un produit ligne block
                else {
                    cout << "la" << endl;
                    const Matrix<double> *ktemp = H->get_dense_block_data();
                    Matrix<double> kkk          = *ktemp;
                    for (int l = 0; l < kkk.nb_rows(); ++l) {
                        const vector<double> xk = kkk.get_row(l);
                        vector<double> yk(K->get_source_cluster().get_size(), 0);
                        produitP(H, xk, yk, 'T');
                        cout << "prod ok" << endl;
                        for (int ii = 0; ii < K->get_source_cluster().get_size(); ++ii) {
                            Mtemp(l, ii) += yk[ii];
                        }
                    }
                    cout << "assemblage ok" << endl;
                    // val = val + temp;
                }
            }
        }
        //Matrix<double>* val = &Mtemp;
        const Matrix<double>* val = Lk->get_dense_block_data();
        val = &Mtemp;
        }
         */
    }
}


int main(int argc, char *argv[]) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Check the number of parameters
    if (argc < 1) {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] << " outputpath" << endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }

    std::string outputpath = argv[1];

    // Htool parameters
    double epsilon = 0.001;
    double eta     = 100;

    // n² points on a regular grid in a square
    int n    = std::sqrt(4761);
    int size = n * n;

    // p1: points in a square in the plane z=z1
    double z = 1;
    vector<double> p(3 * size);
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            p[3 * (j + k * n) + 0] = j;
            p[3 * (j + k * n) + 1] = k;
            p[3 * (j + k * n) + 2] = z;
        }
    }

    // Hmatrix
    MyMatrix A(3, size, size, p, p);
    std::vector<double> x(size, 1), result(size, 0);
    std::shared_ptr<Cluster<PCA<SplittingTypes::RegularSplitting>>> t = make_shared<Cluster<PCA<SplittingTypes::RegularSplitting>>>(3);
    t->build(size, p.data(), 2);
    HMatrix<double> HA(t, t, epsilon, eta, 'N', 'N');
    HA.build(A, p.data());
    result = HA * x;

    // Output
    HA.print_infos();
    HA.save_plot(outputpath + "/smallest_example_plot");
    HA.get_target_cluster()->save_geometry(p.data(), outputpath + "/smallest_example_cluster", {1, 2, 3});
    std::cout << outputpath + "/smallest_example_plot" << std::endl;
    std::cout << Frobenius_absolute_error(HA, A) / A.norm() << std::endl;
    std::cout << norm2(A * x - result) / norm2(A * x) << std::endl;
    cout<<"test"<<endl;

    // On rajoute nos test.


    // On veut modifier les blocs de L =HK

    // On copie L ( j'arrive pas a l'initialiser a 0
    HMatrix<double> Res(t, t, epsilon, eta, 'S', 'U');
    Res.build(A, p.data());
    vector<Block<double>*> Data = HA.get_ComputedBlock();
    vector<Block<double>*> DataRes = Res.get_ComputedBlock();
    // On va classer les low rank et pas low rank


//     vector<Block<complex<double>>> ll,fl,llR,flR;

//     for(auto it = Data.begin(); it!= Data.end(); ++it){
//       int ofs_t =  it->get_target_cluster.get_offset(); int ofs_s  = it->get_source_cluster.get_offset();
//       int sz_t =  it->get_target_cluster.get_size(); int sz_s  = it->get_source_size.get();
//       if( 
// fonction restrict
    // int eps0,eps1,eps2;
    // eps0=0;eps1=0;eps2=0;
    // int nr = 0 ;int nf = 0;
    // for(auto it = Data.begin(); it!= Data.end();++it){
    //   int temp = 0;
    //   Block<double>* bb = *it;
    //   int sz_t = bb->get_target_cluster().get_size();int sz_s = bb->get_source_cluster().get_size();
    //   int of_t = bb->get_target_cluster().get_offset(); int of_s = bb->get_source_cluster().get_offset();
    //   vector<double> x1(sz_s,1); vector<double> x2(sz_s,1);
    //   vector<double> y1(sz_t,0); vector<double> y2(sz_t,0);
    //   bb->get_block_data()->add_mvprod_row_major(x1.data(),y1.data(),1,'N','N');
    //   Matrix<double> m (sz_t,sz_s);
    //   for(int k =0; k<sz_t;++k){
    // 	for(int l =0; l<sz_s;++l){
    // 	  m(k,l)=A.get_coef(k+of_t,l+of_s);} }
    //   vector<double> yres = m*x1;
    //   for(int k=0; k<sz_t;++k){
    // 	temp+=((yres[k]-y1[k])*(yres[k]-y1[k]));
    //   }
    //   bool test = (bb->get_low_rank_block_data()==nullptr);
    //   cout<< yres[1]<<","<<test<<","<< sqrt(temp)/sz_t<< endl;
    //   eps0+= temp;}
    // cout<< sqrt(eps0)/4761<<endl;
    vector<double> xx(4761,1); vector<double> yy(4761,0); vector<double> ref(4761,0);
    ref = A*xx;
    Block<double>* bb0 = Data[0]; Block<double>* Rt = bb0->get_root();
    vector<double>* py = &yy;
    int rr = 0;
    Produit(Rt,xx,py,&rr);
    vector<double> teest = *py;
     int nr = 0;
     //for(int k =0; k< 4761; ++k){
       //nr+= (ref[k]-teest[k])*(ref[k]-teest[k]);}
     //cout<< nr<< ','<< sqrt(nr)/4761<<endl;
     cout<<norm2(ref-teest)/norm2(ref)<< endl;
     cout<< rr<< endl;
/*
     //Test pour restrict
     cout<< "Test Restrict"<< endl;
     pair<vector<Matrix<double>>,vector<int>> sr; vector<Block<double>*> sh;
     sh.push_back(Rt); sh.push_back(Rt);
     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> S;
     S.first = sr; S.second=sh;
     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> Stest =  Restrict0(S,Rt->get_target_cluster(),Rt->get_source_cluster());
     pair<vector<Matrix<double>>,vector<int>> srtest = Stest.first; vector<Block<double>*> ssh = Stest.second;
     cout<< (srtest.first).size()<<','<<(srtest.second).size()<<","<<endl;
     cout<< ssh.size()<<','<< ssh[0]->get_target_cluster().get_offset()<<','<<sh[0]->get_target_cluster().get_size()<<endl;
     cout<<".............."<<endl;

     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> Stest1 =  Restrict0(S,Rt->get_target_cluster().get_son(0),Rt->get_source_cluster().get_son(0));
     pair<vector<Matrix<double>>,vector<int>> srtest1 = Stest1.first; vector<Block<double>*> ssh1 = Stest1.second;
     cout<< (srtest1.first).size()<<','<<(srtest1.second).size()<<","<< endl;
     cout<< ssh1.size()<<','<< ssh1[0]->get_target_cluster().get_offset()<<','<<ssh1[0]->get_target_cluster().get_size()<<endl;
     cout<<".............."<<endl;

     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> Stest2 =  Restrict0(S,Rt->get_son(0).get_son(0).get_target_cluster(),Rt->get_son(0).get_son(2).get_source_cluster());
     pair<vector<Matrix<double>>,vector<int>> srtest2 = Stest2.first; vector<Block<double>*> ssh2 = Stest2.second;
     cout<< (srtest2.first).size()<<','<<(srtest2.second).size()<<endl;
     cout<< ssh2.size()<<','<< ssh2[0]->get_target_cluster().get_offset()<<','<<ssh2[0]->get_target_cluster().get_size()<<endl;
     cout<<Rt->get_son(0).get_son(2).get_source_cluster().get_size()<< endl;
     cout<<Rt->get_son(0).get_son(2).get_source_cluster().get_offset()<< endl;
     cout<<".............."<<endl;

*/
     cout<< "Test Restrict0"<< endl;
     pair<vector<Matrix<double>>,vector<int>> sr; vector<Block<double>*> sh;
     Block<double>& B= Rt->get_son(0);Block<double>* BB = &B;
     sh.push_back(BB); sh.push_back(BB);
     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> S;
     S.first = sr; S.second=sh;
     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> Stest =  Restrict0(S,BB->get_target_cluster().get_son(0),BB->get_source_cluster().get_son(1));
     pair<vector<Matrix<double>>,vector<int>> srtest = Stest.first; vector<Block<double>*> ssh = Stest.second;
     cout<< (srtest.first).size()<<','<<(srtest.second).size()<<endl;
     cout<< ssh.size()<<endl;
     cout<< ssh[0]->get_target_cluster().get_offset()<<','<<ssh[0]->get_target_cluster().get_size()<<endl;
     for ( int k =0; k< 2; ++k){
         cout<<"++++++++++++++++"<<endl;
         cout<< ssh[2*k]->get_target_cluster().get_offset()<<','<<ssh[2*k]->get_target_cluster().get_size()<<endl;
         cout<< ssh[2*k]->get_source_cluster().get_offset()<<','<<ssh[2*k]->get_source_cluster().get_size()<<endl;
         cout<<"------------------"<< endl;
             cout<< ssh[2*k+1]->get_target_cluster().get_offset()<<','<<ssh[2*k+1]->get_target_cluster().get_size()<<endl;
         cout<< ssh[2*k+1]->get_source_cluster().get_offset()<<','<<ssh[2*k+1]->get_source_cluster().get_size()<<endl;


     }
     cout<<".............."<<endl;
    // test avec les valeurs des blocs
    Matrix<double> M0(Rt->get_target_cluster().get_size(),Rt->get_source_cluster().get_size());
    for(int k =0; k < Rt->get_target_cluster().get_size(); ++k){
        for (int l = 0; l< Rt->get_source_cluster().get_size(); ++l){
            M0(k,l)=A.get_coef(k,l);
        }
    }
     cout<<"test blocs"<< endl;
     vector<int> rep;
     for (int k =0; k< Data.size(); ++k){
         Block<double>* bb = Data[k];
         if((bb->get_target_cluster().get_offset()==bb->get_source_cluster().get_offset())and (bb->get_target_cluster().get_size()==bb->get_source_cluster().get_size())){
             rep.push_back(k);
         }
     }
     int cpt = 0; int cps = 0;
       //on en predn un qui peut etre mis au carré juste pour testé
     for (int i = 0; i < rep.size() ; ++i) {
         int kk            = rep[i];
         Block<double> *bb = Data[kk];

         int i1                                = bb->get_target_cluster().get_size();
         int i2                                = bb->get_source_cluster().get_size();
         int o1                                = bb->get_target_cluster().get_offset();
         int o2                                = bb->get_source_cluster().get_offset();
         const VirtualBlockData<double> *vtest = bb->get_block_data();
         // double test = coef(0,0,i1,i2,bb0->get_block_data());
         cout << coef(1, 1, i1, i2, vtest) << "           et          " << M0(o1+1, o2+1) << endl;
         cout<< o1<<','<<o2<<endl;
         Matrix<double> mm(i1, i2);

         for (int k = 0; k < i1; ++k) {
             //cout<<"-----------------"<<endl;
             for (int l = 0; l < i2; ++l) {
                 mm(k, l) = M0(k+o1,l+o2);
                 //if (coef(k,l, i1, i2, vtest)  == A.get_coef(o1+k, o2+l) ){ cpt+=1;}
                 //else{cps+=1;}

             }
         }
         //vector<double> yy = Col(0,i1,i2,vtest);cout<<yy[0]<<','<<yy[1]<<endl;
         vector<double> yres(Rt->get_target_cluster().get_size(),0);
         vector<double> x0(Rt->get_source_cluster().get_size(),1);
         vector<double> yyres(bb->get_target_cluster().get_size(),0);
         vector<double> xtest(i2,1);
         yyres = mm*x;

         vector<double>* Y = &yres;
         int rr0 = 0;
         Produit(bb,x0,Y,&rr);
         vector<double> test0= *Y;
         vector<double> ytemp;
         for(int k0=0 ; k0< i2; ++k0){
             ytemp.push_back(test0[bb->get_source_cluster().get_offset()+k0]);
         }
         cout<< "++++++++"<< endl;
         //cout<<norm2(yy-mm.get_col(0))<<endl;
         cout<<norm2(ytemp-yyres)/norm2(yyres)<< endl;

         //cout<<cpt<<','<<cps<< endl;
         //Matrix<double> test = get_mat(i1, i2, vtest, vtest);
         //Matrix<double> ref  = Mprod(mm,mm);
         //cout << normFrob(test-ref)<<"!!!!"<< normFrob(test-ref)/ normFrob(ref)<<endl;
     }
     //cout<<cpt<<','<<cps<< ','<<cpt+cps<<endl;
    // cout<< Rt->get_source_cluster().get_son(0).get_size()<<","<< Rt->get_source_cluster().get_son(0).get_offset()<< endl;
    //cout<< Rt->get_target_cluster().get_size()*Rt->get_source_cluster().get_size()<<endl;
     //cout<<Rt->nb_sons()<<endl;
     //cout<< Rt->get_son(0).nb_sons()<< endl;
     //Rt->get_data

/*
         Block mult(block K){

         Block L(this->admissibility_condition,this->t,K.get_source_cluster(),this,K);
         std::stack<Block*> stack_block;
         std::stack<SumExpression*> stack_sum_expression;
         stack_block.push(L);
         while (!stack_block.empty()){
             current_block = stack_block.top();
             current_sum_expression = stack_sum_expression.top();

             if (!(current_block->get_target_cluster()->IsLeaf()) && !(current_block->get_source_cluster()->IsLeaf())){

                 for (int i=0;i<current_block->get_target_cluster()->nb_sons();i++){
                     for (int j=0;j<current_block->get_source_cluster()->nb_sons();j++){

                         current_block->sons.push_back(std::unique_ptr<Block>(new Block(this->admissibility_condition, current_block->get_target_cluster()get_son(i), current_block->get_source_cluster()->get_son(j), &L, L.tasks, L.local_tasks)));

                         current_block->sons.back()->sum_expression =  current_block.sum_expression.restrict(i,j);

                         // SumExpression
                     }
                 }

                 for (int p=0;p<current_block->sons.size();p++){
                     stack_block.push(current_block->sons[p]);
                 }
             }
             else{
                 bool Admissible= this->admissibility_condition(this->t,K.get_source_cluster());
                 if(Admissible){
                     lrmat_block_data = new DenseBlockData<T>()
                                            current_sum_expression.truncate(lrmat_block_data);

                     block_data = std::unique_ptr<VirtualBlockData<T>>(low_rank_block_data);
                 }
                 else{
                     // Possibilité 1
                     dense_block_data = new DenseBlockData<T>(current_sum_expression);

                     // Possibilité 2
                     dense_block_data = new DenseBlockData<T>();
                     current_sum_expression.evaluate(dense_block_data);

                     //
                     current_block->block_data = std::unique_ptr<VirtualBlockData<T>>(dense_block_data);

                 }

             }

         }

         return L;
     }
*/
//Test Produit

     cout <<"Test Produit"<<endl;
     vector<double> x0(4761,2), yy0(4761,0), ref0(4761,0),xx_perm(4761,0),yy_perm(4761,0);
     global_to_cluster(&Rt->get_source_cluster(),xx.data(),xx_perm.data());
     cluster_to_global(&Rt->get_target_cluster(),yy0.data(),yy_perm.data());
     produitP(Rt,xx_perm,yy_perm,'T');
     cluster_to_global(&Rt->get_target_cluster(),yy_perm.data(),yy0.data());
     //vector<double> yref0 = M0*x0;

     vector<double> yref0;
     for (int l =0 ; l < Rt->get_source_cluster().get_size(); ++l){
        double tt = 0;
        for (int k =0; k< Rt->get_target_cluster().get_size();++k){
            tt+=  M0(k,l)*x[k];
     }
     yref0.push_back(tt);}
     cout << yref0.size()<< endl;
     cout<< norm2(yy0-yref0)/norm2(yref0)<<endl;
     //Produit marche pour les deux sens !!!!
     /*
     for(int k =0; k < Data.size(); ++k){
        Block<double>* bb0 = Data[k]; Block<double>* Rt = bb0->get_root();
        Produit(Rt,xx,py,&rr);
     //for(int k =0; k< 4761; ++k){
     //nr+= (ref[k]-teest[k])*(ref[k]-teest[k]);}
     //cout<< nr<< ','<< sqrt(nr)/4761<<endl;
     cout<<norm2(ref-teest)/norm2(ref)<< endl;
     cout<< rr<< endl;
     */
     //TEST SUM EXPRESSION
     cout<<"___________________"<< endl;
     SumExpression Sumex(&Rt->get_son(0),&Rt->get_son(0));
     cout << Sumex.get_sh().size()<<','<< Sumex.get_sr().size()<<endl;


    SumExpression St = Sumex.Restrict(Rt->get_son(0).get_target_cluster().get_son(0),Rt->get_son(0).get_source_cluster().get_son(1));
    cout << St.get_sh().size()<<','<< St.get_sr().size()<<endl;
    vector<Block<double>*> bbb = St.get_sh();
    SumExpression Stt = Sumex.Restrict(Rt->get_son(0).get_target_cluster().get_son(0),Rt->get_son(0).get_source_cluster().get_son(1));
    cout << Stt.get_sh().size()<<','<< Stt.get_sr().size()<<endl;
    vector<Block<double>*> bbbt = Stt.get_sh();
    for ( int k =0 ; k < Stt.get_sh().size();++k){
        Block<double>* bt = bbbt[k];
        Block<double>* b = bbb[k];

        cout<<"///////"<<endl;
        cout<< bt->get_target_cluster().get_offset()<<','<<bt->get_target_cluster().get_size()<<endl;
        cout<<bt->get_source_cluster().get_offset()<<','<<bt->get_source_cluster().get_size()<<endl;
        cout<<"=========="<<endl;
        cout<< b->get_target_cluster().get_offset()<<','<<b->get_target_cluster().get_size()<<endl;
        cout<<b->get_source_cluster().get_offset()<<','<<b->get_source_cluster().get_size()<<endl;



    }



cout<<"teeeeeeeeeeeeeeeeest"<<endl;
cout <<(bb0->get_dense_block_data()==nullptr) <<','<<(bb0->get_low_rank_block_data() == nullptr)<<endl;
const Matrix<double>* mmm = bb0->get_dense_block_data();
Matrix<double> m1 = *mmm;
cout<<m1(1,1)<<','<<m1.nb_cols()<<','<<m1.nb_rows()<<endl;

cout<<"test produit hmat"<<endl;

//Matrix<double> m0 = Mprod(M0,M0);
Mycondition *my;
Block<double> root ( my, Rt->get_target_cluster(), Rt->get_source_cluster());
SumExpression Su(Rt,Rt);
//test pour build un son
cout<<root.nb_sons()<<endl;
//vector<unique_ptr<Block<double>>> pt0 = root.get_Son();cout<<pt0.size()<<endl;
root.build_son(Rt->get_target_cluster(),Rt->get_source_cluster());
cout<< root.nb_sons()<<endl;
//auto sons = root.Get_Son();
Hmult(&root,&root,Su);
//Block<double>* Res0 = &root;

/*
vector<double> xtest(4761,2), ytest(4761,0), test(4761,0),xx_test(4761,0),yy_test(4761,0);
global_to_cluster(&Rt->get_source_cluster(),xtest.data(),xx_test.data());
cluster_to_global(&Rt->get_target_cluster(),ytest.data(),yy_test.data());
produitP(Res0,xx_test,yy_test,'N');
cluster_to_global(&Rt->get_target_cluster(),yy_test.data(),ytest.data());
cout<< norm2(M0*(M0*xtest)-yy_test)/norm2(M0*(M0*xtest))<<endl;




*/

MPI_Finalize();
}
