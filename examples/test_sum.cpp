#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <complex>
#include <vector>
#include <htool/htool.hpp>
#include <typeinfo>
using namespace std;
using namespace htool;


// Condition d'adm 
zf zac
class MyCondition: public VirtualAdmissibilityCondition {
   bool ComputeAdmissibility(const VirtualCluster &target, const VirtualCluster &source, double eta) const override {
        bool admissible = 2 * std::min(target.get_rad(), source.get_rad()) < eta * std::max((norm2(target.get_ctr() - source.get_ctr()) - target.get_rad() - source.get_rad()), 0.);
        return admissible;
    }
};


// Produit hmat/bloc vecteur

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


//classe Sum-expr
 class SumExpression {
  private :
   // liste de matrice de rang faible, d'où elles sont sutiué, liste de hmat
    vector<Matrix<double>> SR;
    vector<int> off;
    vector<Block<double>*> SH;

  public:
    //Constructeur qui prend deux blocs
    SumExpression(Block<double>* H, Block<double>* K){ vector<Block<double>*> vb; vb.push_back(H); vb .push_back(K); SH = vb;}
    // Constructeur avec tout
    SumExpression(const vector<Matrix<double>> SR0,const vector<int> off0,const vector<Block<double>*> Sh0){this->SR=SR0;this->off=off0;this->SH=Sh0;}
   //getters
    vector<Matrix<double>> get_sr(){return this->SR;}
    vector<int> get_off(){return this->off;}
    vector<Block<double>*> get_sh(){return this->SH;}

 };


//classe mymatrix pour nos tests
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
};


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

pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict(pair<pair < vector<Matrix<double>>,vector<int>>, vector<Block<double>*>>S,const VirtualCluster& t, const VirtualCluster& s){
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


int main(int argc, char *argv[]) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Check the number of parameters
    if (argc < 1) {
        cerr << "Usage: " << argv[0] << " outputpath" << endl;
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

    Block<double>* bb0 = Data[0]; Block<double>* Rt = bb0->get_root();

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

MPI_Finalize();
}
