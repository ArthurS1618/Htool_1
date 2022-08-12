#include <htool/htool.hpp>


using namespace std;
using namespace htool;

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








    MPI_Finalize();
}
