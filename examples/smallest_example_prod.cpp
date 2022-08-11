#include <htool/htool.hpp>


using namespace std;
using namespace htool;


class MyCondition: public VirtualAdmissibilityCondition {
   bool ComputeAdmissibility(const VirtualCluster &target, const VirtualCluster &source, double eta) const override {
        bool admissible = 2 * std::min(target.get_rad(), source.get_rad()) < eta * std::max((norm2(target.get_ctr() - source.get_ctr()) - target.get_rad() - source.get_rad()), 0.);
        return admissible;
    }
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



};
//   //////////////////////////////////////////////////
  // std::vector<Block<double>> restrict(std::pair<std::vector<Block<double>>,std::vector<pair<Block<double>,Block<double>>>> S, const VirtualCluster& t, const VirtualCluster& s){
  // std::vector<Block<double>> Sr = S.first;
  // //il faut faire la restriction des low rank a t,s
  // std::vector<Block<double>> Sr0;    int offset_t = t.get_offset(); int offset_s = s.get_offset(); int size_t = t.get_size(); int size_s = s.get_size();
  // for (auto it= Sr.begin(); it!= Sr.end(); ++it){
  //   // Block<double> &lr = *it;
  //   int repere = 0;
  //   // // ____________________//
  //   // //Pas sure que cette méthode marche mais sinon il faut mettre la matrice en full, faire sa restriction et refaire ACA;
  //   // // on peut aussi juste faire une boucle for k < lr.get_nb_son()
  //   // //____________________//
  //   // // On récupère la réstriction a tau sigma en vérifiant sur chaques blocs des get_son() si leurs clusters sont les mêmes
  //   // while(!(&lr.get_son(repere) ==nullptr)){
  //   //   Block<double> &lr_son = lr.get_son(repere);
  //   //   int offset_t0=lr_son.get_target_cluster().get_offset(); int offset_s0= lr_son.get_source_cluster().get_offset();
  //   //   int size_t0= lr_son.get_target_cluster().get_size(); int size_s0 = lr_son.get_source_cluster().get_size();
  //   //   // on vérifie que les cluster sont les mêmes <=> même sze et offset
  //   //   if( ((offset_t == offset_t0) and (offset_s == offset_s0))  and ((size_t0 == size_t) and (size_s0 == size_s))){
  //   // 	Sr0.push_back(lr_son);}
  //     repere +=1;
  //   }
  
  // //on initialise Sh0 a 0
  // std::vector<pair<Block<double>,Block<double>>> Sh = S.second;
  // std::vector<pair<Block<double>,Block<double>>> Sh0;
  // //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
  // for (auto it= Sh.begin(); it!= Sh.end(); ++it){
  //   pair<Block<double>,Block<double>> HK = *it; Block<double> H=HK.first; Block<double> K = HK.second;
  //   Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
  //   int repere1 = 0;
  //   //on parcours les enfants de rho=r
  //   while(!(r.get_son_ptr(repere1) == nullptr)){
  //     //on cherche les restrictions de H et K parmi leur fils comme pour Sr
  //     int repereH =0;int repereK=0;Block<double> H0;Block<double> K0;
  //     int offset_r0 = r.get_son_ptr(repere).get_offset();
  //     int size_r0= r.get_son_ptr(repere).get_size();
  //     //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
  //     while(!(&H.get_son(repereH) == nullptr)){
  // 	Block<double> &Htemp = H.get_son(repereH);
  // 	int offset_t0 = Htemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Htemp.get_source_cluster().get_offset();
  // 	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
  // 	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
  // 	  H0 = Htemp;}
  // 	repereH +=1;}
  //     // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
  //     while(!(&K.get_son(repereK) == nullptr)){
  // 	Block<double> &Ktemp = K.get_son(repereK);
  // 	int offset_t0 = Ktemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Ktemp.get_source_cluster().get_offset();
  // 	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
  // 	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
  // 	  K0 = Ktemp;}
  // 	repereK +=1;}
  //     //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
  //     // Si c'est le cas le produit est forcément low rank
  // 	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
  // 	//On veut faire ABt = H0*K0------>push_back->Sr0
  // 	// trois possibilité: lr*lr, lr*hmat, hmat*lr
  // 	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
  // 	// Dans les trois cas on a un lr en résultat
  // 	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
  // 	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
  // 	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
  // 	 LowRankMatrix<double>& hr0 = H0.get_low_rank_block_data();
  // 	 LowRankMatrix<double>& kr0 = K0.get_low_rank_block_data();
  // 	 DenseBlockData<double>& hf0 = H0.get_dense_block_data();
  // 	 DenseBlockData<double>& kf0 = K0.get_dense_block_data();
  // 	// //cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
  // 	// if((hf0 == nullptr) and (kf0 == nullptr)){
  // 	//   LowRankMatrix<T> hk0 = h0*k0;}
  // 	// //cas 2 : hmat*lowrank
  // 	// else if( hr0 == nullptr ){
  // 	//   Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
  // 	//   int rank = U.nb_cols();
  // 	//   LowRankMatrix<T> hk0;
  // 	//   Matrix<T>* hk0U = hk0.Get_U();
  // 	//   Matrix<T>* hk0V = res.Get_V();
  // 	//   *hk0V = V;
  // 	//   Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
  // 	//   for ( int k =0; k< rank;++k){
  // 	//     vector<T> Uk = U.get_col(k);
  // 	//     vector<T> temp = H0*Uk;
  // 	//     for ( int l =0; l< U.nb_rows();++l){
  // 	//       mattemp(l,k)= Uk[l];
  // 	//     }
  // 	//   }
  // 	//   *hk0U = mattemp;
  // 	// }
  // 	// //cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
  // 	// // Sinon on code un produit vecteur mat xA = At x
  // 	// else if( kr0 ==nullptr){
	  
  // 	// }
	  
  //     }
  //     else{
  // 	pair<Block<double>,Block<double>> HK0; HK0.first = H0;HK0.second = K0;
  // 	Sh0.push_back(HK0);
  //     }
  //   }
  // return Sr0;
  //}
  // }//   //////////////////////////////////////////////////




  
//   std::pair<std::vector<Block<double>>,std::vector<pair<Block<double>,Block<double>>>> restrict(std::pair<std::vector<Block<double>>,std::vector<pair<Block<double>,Block<double>>>> S, const VirtualCluster& t, const VirtualCluster& s){
//   std::vector<Block<double>> Sr = S.first;
//   //il faut faire la restriction des low rank a t,s
//   std::vector<Block<double>> Sr0;    int offset_t = t.get_offset(); int offset_s = s.get_offset(); int size_t = t.get_size(); int size_s = s.get_size();
//   for (auto it= Sr.begin(); it!= Sr.end(); ++it){
//     Block<double> &lr = *it;
//     int repere = 0;
//     // ____________________//
//     //Pas sure que cette méthode marche mais sinon il faut mettre la matrice en full, faire sa restriction et refaire ACA;
//     // on peut aussi juste faire une boucle for k < lr.get_nb_son()
//     //____________________//
//     // On récupère la réstriction a tau sigma en vérifiant sur chaques blocs des get_son() si leurs clusters sont les mêmes
//     while(!(&lr.get_son(repere) ==nullptr)){
//       Block<double> &lr_son = lr.get_son(repere);
//       int offset_t0=lr_son.get_target_cluster().get_offset(); int offset_s0= lr_son.get_source_cluster().get_offset();
//       int size_t0= lr_son.get_target_cluster().get_size(); int size_s0 = lr_son.get_source_cluster().get_size();
//       // on vérifie que les cluster sont les mêmes <=> même sze et offset
//       if( ((offset_t == offset_t0) and (offset_s == offset_s0))  and ((size_t0 == size_t) and (size_s0 == size_s))){
// 	Sr0.push_back(lr_son);}
//       repere +=1;
//     }
//   }
//   // //on initialise Sh0 a 0
//   std::vector<pair<Block<double>,Block<double>>> Sh = S.second;
//   std::vector<pair<Block<double>,Block<double>>> Sh0;
//   // //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
//   // for (auto it= Sh.begin(); it!= Sh.end(); ++it){
//   //   pair<Block<double>,Block<double>> HK = *it; Block<double> H=HK.first; Block<double> K = HK.second;
//   //   Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
//   //   int repere1 = 0;
//   //   //on parcours les enfants de rho=r
//   //   while(!(r.get_son_ptr(repere1) == nullptr)){
//   //     //on cherche les restrictions de H et K parmi leur fils comme pour Sr
//   //     int repereH =0;int repereK=0;Block<double> H0;Block<double> K0;
//   //     int offset_r0 = r.get_son_ptr(repere).get_offset();
//   //     int size_r0= r.get_son_ptr(repere).get_size();
//   //     //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
//   //     while(!(&H.get_son(repereH) == nullptr)){
//   // 	Block<double> &Htemp = H.get_son(repereH);
//   // 	int offset_t0 = Htemp.get_target_cluster().get_offset();
//   // 	int offset_s0 = Htemp.get_source_cluster().get_offset();
//   // 	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
//   // 	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
//   // 	  H0 = Htemp;}
//   // 	repereH +=1;}
//   //     // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
//   //     while(!(&K.get_son(repereK) == nullptr)){
//   // 	Block<double> &Ktemp = K.get_son(repereK);
//   // 	int offset_t0 = Ktemp.get_target_cluster().get_offset();
//   // 	int offset_s0 = Ktemp.get_source_cluster().get_offset();
//   // 	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
//   // 	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
//   // 	  K0 = Ktemp;}
//   // 	repereK +=1;}
//   //     //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
//   //     // Si c'est le cas le produit est forcément low rank
//   // 	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
//   // 	//On veut faire ABt = H0*K0------>push_back->Sr0
//   // 	// trois possibilité: lr*lr, lr*hmat, hmat*lr
//   // 	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
//   // 	// Dans les trois cas on a un lr en résultat
//   // 	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
//   // 	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
//   // 	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
//   // 	 LowRankMatrix<double>& hr0 = H0.get_low_rank_block_data();
//   // 	 LowRankMatrix<double>& kr0 = K0.get_low_rank_block_data();
//   // 	 DenseBlockData<double>& hf0 = H0.get_dense_block_data();
//   // 	 DenseBlockData<double>& kf0 = K0.get_dense_block_data();
//   // 	// //cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
//   // 	// if((hf0 == nullptr) and (kf0 == nullptr)){
//   // 	//   LowRankMatrix<T> hk0 = h0*k0;}
//   // 	// //cas 2 : hmat*lowrank
//   // 	// else if( hr0 == nullptr ){
//   // 	//   Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
//   // 	//   int rank = U.nb_cols();
//   // 	//   LowRankMatrix<T> hk0;
//   // 	//   Matrix<T>* hk0U = hk0.Get_U();
//   // 	//   Matrix<T>* hk0V = res.Get_V();
//   // 	//   *hk0V = V;
//   // 	//   Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
//   // 	//   for ( int k =0; k< rank;++k){
//   // 	//     vector<T> Uk = U.get_col(k);
//   // 	//     vector<T> temp = H0*Uk;
//   // 	//     for ( int l =0; l< U.nb_rows(//   // 	//       mattemp(l,k)= Uk[l];
//   // 	//     }
//   // 	//   }
//   // 	//   *hk0U = mattemp;
//   // 	// }
//   // 	// //cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
//   // 	// // Sinon on code un produit vecteur mat xA = At x
//   // 	// else if( kr0 ==nullptr){
	  
//   // 	// }
	  
//   //     }
//   //     else{
//   // 	pair<Block<double>,Block<double>> HK0; HK0.first = H0;HK0.second = K0;
//   // 	Sh0.push_back(HK0);
//   //     }
//   //   }
//   // }
//   std::pair<std::vector<Block<double>>,std::vector<pair<Block<double>,Block<double>>>> S0;
//   S0.first = Sr0; S0.second = Sh0;
//   return S0;
// }  

//   ////////////////////////////////////



/////////////////// Ca ca marche pour les low rankr

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
	  sk.push_back(Uk);
	  sk.push_back(Vk);
	  ofk.push_back(ofu+repere_i);ofk.push_back(ofv+repere_j);}
  }
    SR.first = sr ; SR.second = of;

  return SR;}

// pair<pair < vector<Matrix<double>> ,vector<int> >,vector<pair< const VirtualCluster&, const VirtualCluster&>>> Restrict (pair<pair < vector<Matrix<double>> ,vector<int> >,vector<pair< const VirtualCluster&,const  VirtualCluster&>>>S, const VirtualCluster& t, const VirtualCluster &s){
// //   //pair<vector<Matrix<double>>,vector<int>> SH = S.first;

// return S;}
//pair<pair < vector<Matrix<double>> ,vector<int> >,vector<VirtualCluster*>> Restrict (pair<pair < vector<Matrix<double>> ,vector<int> >,vector<VirtualCluster*>>S, const VirtualCluster& t, const VirtualCluster &s){
//   //pair<vector<Matrix<double>>,vector<int>> SH = S.first;

//return S;}
////////////////
  
double NORM(vector<double>x,vector<double>y){
  double nr= 0;
  for(int k =0; k < x.size();++k){
    nr+=((x[k]-y[k])*(x[k]-y[k]));
  }
  return sqrt(nr);}


// void mvprod(Block<double>* B,vetcor<double>*y , vector<double>x , vector<Block<double>*> Tb){
//   int of_t = B->get_target_cluster().get_offset(); int sz_t = B->get_target_cluster().get_size();
//   int of_s = B->get_source_cluster().get_offset();int sz_s = B->get_source_clustr().get_size();
//   for (int k =0; k< Tb.size(); ++k){
// 	Block<double>* bk = Tb[k];
// 	 int ofk_t = bk->get_target_cluster().getsize(); int szk_t = bk->get_target_cluster().get_size();
// 	 int ofk_s = bk->get_source_cluster().get_offset();int szk_s = bk->get_source_clustr().get_size();
// 	 if ((of_t==ofk_t) and ( of_s==ofk_s) and (sz_t==szk_t) and (sz_s==szk_s)){

// 	      vector<double> xk;
// 	      for( int i = 0; i<sz_s; ++i){
// 		xk.push_back(x[B->get_source_cluster().get_offset()+i]);}
	      
// 	       bk->get_block_data()->add_mvprod_row_major(y + (ofk_t -z) * mu, temp.data() + offset_j * mu, mu, 'N', 'N');;
// 	      vector<double> yy = temp1;
// 	      vector<double> py = *y;
// 	      vector<double> test(py.size(),0);

//       *y=test;   }  
	   
//       }}
//     else{}}

void Produit(Block<double>*B,vector<double>* y, vector<double> x){
  if (B->nb_sons() == 0){
    if( (B->get_low_rank_block_data() == nullptr) and( B->get_dense_block_data()==nullptr) and(B->get_block_data() ==nullptr) ){
      //cout<< "cheeeeeeeeeeeeeeeeeeeeeelou"<<endl;
      //cout<<B->nb_sons()<<endl;
      //cout<<B->IsAdmissible()<<endl;
      // if(B->get_local_diagonal_block()==nullptr){
      // 	nn=nn+1;
      //   cout<<"!" <<nn<<endl;
      //   //cout<<"la je comprend pas"<<endl;
      // 	//cout<<"..........."<<endl;
      // 	//cout<< B->get_rank_of()<<endl;
      //cout<<B->get_target_cluster().get_offset()<<","<<B->get_target_cluster().get_size()<<endl;
      // 	//cout<<B->get_source_cluster().get_offset()<<","<<B->get_source_cluster().get_size()<<endl;
      // }
      cout<< B->get_target_cluster().get_size()<<","<< B->get_target_cluster().get_offset()<<endl;
      cout<< B->get_source_cluster().get_size()<<","<< B->get_source_cluster().get_offset()<<endl;
      //cout<<1111111111111<< endl;
      }
    
    else if( B->IsAdmissible() ){
      if (!(B->get_low_rank_block_data() ==nullptr)){
      Matrix<double> U = B->get_low_rank_block_data()->Get_U();
      Matrix<double> V = B->get_low_rank_block_data()->Get_V();
      vector<double> xk;
       for( int i = 0; i< V.nb_rows(); ++i){
      	xk.push_back(x[B->get_source_cluster().get_offset()+i]);}
       //vector<double> yy = U*(V*xk);
      vector<double> temp(V.nb_rows());
      for (int k =0; k < V.nb_rows();++k){
	double val = 0;
	for(int l =0; l<V.nb_cols(); ++l){
	  val+=V(k,l)*x[l];}
	temp[k]=val;}
      vector<double> temp1(U.nb_rows());
      for (int k =0; k < U.nb_rows();++k){
	double val = 0;
	for(int l =0; l<U.nb_cols(); ++l){
	  val+=U(k,l)*temp[l];}
	temp1[k]=val;}
      vector<double> yy = temp1;
      vector<double> py = *y;
      vector<double> test(py.size(),0);
      for (int rep = 0; rep < py.size() ; ++rep){
      	if( ((B->get_target_cluster().get_offset()-rep+1)>0) and ((rep-B->get_target_cluster().get_offset())<yy.size())){
      	double yk = py[rep];
        double res = yk+ yy[rep-B->get_target_cluster().get_offset()];
	
	test[rep]=res;
      	// py[rep] = res;
      	  }
      	  else{test[rep] = py[rep];}}
      *y=test;  } }
      // for(int k =0; k< B->get_target_cluster().get_size();++k){
      //	double yk = py[k+B->get_target_cluster().get_offset()];
      // 	double ytemp = yk; double res = ytemp + yy[k];
      // 	py[k+B->get_target_cluster().get_offset()] = res;
      // }
      // y = &py;
      // }}
    else{
      if( !(B->get_dense_block_data() ==nullptr)){
      const Matrix<double>* M = B->get_dense_block_data();
      	Matrix<double> m = *M;
	vector<double>  xk;
        for(int k = 0; k< B->get_source_cluster().get_size(); ++k){
	  xk.push_back(x[B->get_source_cluster().get_offset()+k]);}
  	vector<double> yy = m*xk;
	vector<double> py = *y;
	vector<double> test(py.size(),0);
	for (int rep = 0; rep < y->size() ; ++rep){
	  if( ((rep-B->get_target_cluster().get_offset())< yy.size()) and ((rep-B->get_target_cluster().get_offset()+1)>0)){
	double yk = py[rep];
	double res = yk+ yy[rep-B->get_target_cluster().get_offset()];
	test[rep]=res;
	 py[rep] = res;
	 //cout<<py[rep]<<','<<res<<","<<test[rep]<<endl;
	  }
	  else{test[rep] = py[rep];}}
	*y=test;}
	 //y= &py;
	 //vector<double> gy = *y;
	 //cout<<'?'<<gy[3+B->get_target_cluster().get_size()]<<','<<py[3+B->get_target_cluster().get_size()]<<','<<res<<endl;
 } } 
  else {
    for (int r = 0; r <B->nb_sons();++r){
      Block<double>& Br = B->get_son(r);
      Produit(&Br,y,x);}
  }
}


	 

//Fonction multiplication Block vecteur
vector<double> Prod(Block<double>* B, vector<double> x,vector<double>y){
  int of_s= B->get_source_cluster().get_offset(); int sz_s = B->get_source_cluster().get_size();
  int of_t= B->get_target_cluster().get_offset(); int sz_t = B->get_target_cluster().get_size();
  if (B->nb_sons()==0 ){
    if( B->IsAdmissible()){//Bloc low rank
    // safe guard mais normalement c'est lr
    if(!(B->get_low_rank_block_data() == nullptr)){
      const LowRankMatrix<double>* Lr= B->get_low_rank_block_data();
      vector<double>  xk;
      for(int k = 0; k< sz_s; ++k){
       	xk.push_back(x[of_s+k]);}
      
       Matrix<double> U = B->get_low_rank_block_data()->Get_U(); Matrix<double> V = B->get_low_rank_block_data()->Get_V();
       vector<double> yk = U*(V*xk);
       for(int l =0; l < sz_t; ++l){
	  y[of_t+l]=y[of_t+l]+ yk[l];
	}
    } }
    else{//bloc dense
      // safeguard
      if(!(B->get_dense_block_data()== nullptr)){
	const Matrix<double>* M = B->get_dense_block_data();
	Matrix<double> m = *M;
	vector<double>  xk;
        for(int k = 0; k< sz_s; ++k){
	   xk.push_back(x[of_s+k]);}
	vector<double> yk = m*xk;
	for(int k =0; k < sz_t; ++k){
	  y[of_t+k] =y[of_t+k]+ yk[k];}
      } } }
    else{
      // du coup on a des sons
      for (int k=0; k< B->nb_sons(); ++k){
    	Block<double>& bk = B->get_son(k);
	vector<double> zk (4761,0.);
        vector<double> yk =Prod(&bk,x,zk);
	 for (int kk = 0; kk< bk.get_target_cluster().get_size(); ++kk){
	  y[kk+bk.get_target_cluster().get_offset()]+=yk[kk+bk.get_target_cluster().get_offset()];}
         //y=Prod(&bk,x,zk);
      //else{
      //	Block<double>&bk = B->get_son(0); return Prod(&bk,x,y);}
      }}
  return y;
}



// //Fonction multiplication Block vecteur
// vector<double> Prod(Block<double>* B, vector<double> x,vector<double>y){
//   int of_s= B->get_source_cluster().get_offset(); int sz_s = B->get_source_cluster().get_size();
//   int of_t= B->get_target_cluster().get_offset(); int sz_t = B->get_target_cluster().get_size();
//   if (B->nb_sons()==0 ){
//     if( B->IsAdmissible()){//Bloc low rank
//     // safe guard mais normalement c'est lr
//     if(!(B->get_low_rank_block_data() == nullptr)){
//       const LowRankMatrix<double>* Lr= B->get_low_rank_block_data();
//       vector<double>  xk;
//       for(int k = 0; k< sz_s; ++k){
//        	xk.push_back(x[of_s+k]);}
      
//        Matrix<double> U = B->get_low_rank_block_data()->Get_U(); Matrix<double> V = B->get_low_rank_block_data()->Get_V();
//        vector<double> yk = U*(V*xk);
//        for(int l =0; l < sz_t; ++l){
// 	 // double* tp = y[k];
// 	 // double ttp = *tp+yk[k];
// 	 // tp = &ttp;
//        	y[of_t+l]+= yk[l];
// 	}
//     } }
//     else{//bloc dense
//       // safeguard
//       if(!(B->get_dense_block_data()== nullptr)){
// 	const Matrix<double>* M = B->get_dense_block_data();
// 	Matrix<double> m = *M;
// 	vector<double>  xk;
//         for(int k = 0; k< sz_s; ++k){
// 	   xk.push_back(x[of_s+k]);}
// 	vector<double> yk = m*xk;
// 	for(int k =0; k < sz_t; ++k){
// 	  y[of_t+k] += yk[k];}
//       } } }
//     else{
//       // du coup on a des sons
//       for (int k=0; k< B->nb_sons(); ++k){
//     	Block<double>& bk = B->get_son(k);
// 	//vector<double>yk = Prod(&bk,x,y);
// 	//for (int kk =0; kk<bk.get_target_cluster().get_size();++kk){
// 	// y[bk.get_target_cluster().get_offset()+kk]+=yk[kk];}
// 	// vector<double>z =Prod( &bk,x,y);
// 	// int of = bk.get_target_cluster().get_offset();
// 	// int sz_t = bk.get_target_cluster().get_size();
// 	// for (int l =0; l <sz_t;++l){
// 	//     y[l+bk.get_target_cluster().get_offset()]+=z[l+of];
// 	//  }
// 	y = Prod(&bk,x,y);
//       //else{
//       //	Block<double>&bk = B->get_son(0); return Prod(&bk,x,y);}
//       }}
//   return y;
// }

// vector<double> pp(Matrix<double> M,int , vector<double>x,vector<double> y){
//   if(M.nb_cols()>10){
//     for (int i =0; i<2){
//       for(int j =0; j<2; +j){
// 	Matrix<double> Mi(M.nb_rows()/2,M.nb_cols()/2);
// 	for(int k =0; k < M.nb_rows()/2; ++k){
// 	  for (int l =0; l< M.nb_cols();++l){
// 	    Mi(k,l)=M(i*M.nb_rows()/2,j*M.nb_cols()/2);
// 	  }
// 	}
// 	y = pp(Mi,x,y)
	

int fac(int k ){
  if (k==0){
    return 1;}
  else{return k*fac(k-1);}}
///
pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Restrict(pair < vector<Matrix<double>>,vector<int>>SR0, vector<Block<double>*>SH0,const VirtualCluster& t, const VirtualCluster& s){
  int of_t = t.get_offset(); int sz_t = t.get_size(); int of_s = s.get_offset(); int sz_s = s.get_size();
  pair<vector<Matrix<double>>,vector<int>> Sr = restrict_lr(SR0,t,s);
  pair<vector<Matrix<double>>,vector<int>> SR;
  pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> Res;
  Res.first = SR0;Res.second= SH0;
  int n = SH0.size()/2;
  // On parcours les HK de SH
  for (int k =0; k< n; ++k){
    Block<double>*H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
    //On prend rho=r = H.source= K.target();
    int of_r = K->get_target_cluster().get_offset(); int sz_r = K->get_target_cluster().get_size();
    int nb_r = K->get_target_cluster().get_nb_sons();
    int of_h = H->get_target_cluster().get_offset(); int sz_h =H->get_target_cluster().get_size();
    int of_k = K->get_source_cluster().get_offset(); int sz_k = K->get_source_cluster().get_size();
    //Cas ou son.target=! t
    if ( ((sz_h == sz_t ) and (sz_k== sz_s)) and (( of_h  == of_t) and (of_k == of_s))){
      for (int i = 0; i< nb_r; ++i){
	// On test si l'un des deux est une feuille.
	if( (H->IsAdmissible()) or ( K->IsAdmissible())){
	  
	  //trois cas
	  if( !(H->get_low_rank_block_data() == nullptr) and !( K->get_low_rank_block_data()==nullptr)){
	    // Les deux sont low rank
	    Matrix<double> Uh = H->get_low_rank_block_data()->Get_U(); Matrix<double> Vh = H->get_low_rank_block_data()->Get_V();
	    Matrix<double> Uk = K->get_low_rank_block_data()->Get_U(); Matrix<double> Vk = K->get_low_rank_block_data()->Get_V();
	    Matrix<double> u = Uh; Matrix<double> v= Vh*(Uk*Vk);
	    cout<<"!!!"<<endl;
	    cout<< Sr.first.size()<<endl;
	    Sr.first.push_back(u);Sr.first.push_back(v);Sr.second.push_back(of_t);Sr.second.push_back(of_s);
	    cout<< Sr.first.size()<<endl;
	   	}
	  else if( !(K->get_low_rank_block_data()==nullptr)){
	    Matrix<double> Uk = K->get_low_rank_block_data()->Get_U();Matrix<double> Vk = K->get_low_rank_block_data()->Get_V();
	    int nc = Uk.nb_cols();
	    Matrix<double> (sz_t,nc);
	    for (int r1 = 0; r1 < nc ; ++r1){
	      vector<double> x = Uk.get_col(r1);
	      //La il faut faire la multiplication hmat vecteur mais on a que un bloc !
	      //vector<double> y = 
	    }}
	}
      }
    }
    else{
      cout<<"chelou"<<endl;}

  }
  //for (int k =0; k < n ; endl){
  //Block<double>* H = SH0[2*k]; Block<double>* K = SH0[2*k+1];
    
  //}
  return Res;}




// pair<pair<vector<Matrix<double>>,vector<int>>,vector<pair<VirtualCluster , VirtualCluster>> Restrict(pair<pair<vector<Matrix<double>>,vector<int>>,vector<pair<VirtualCluster,VirtualCluster>>> S, const VirtualCluster &t, const VirtualCluster &s){
//   pair<vector<Matrix<double>>,vector<int>> SR0 , SR;
//   //vector<HMatrix<double>> SH = S.second;
//   SR0 = S.first;
//   SR = restrict_lr(SR0,t,s);
//   //auto tt = S.second;
//   //cout<< typeid(tt).name()<< endl;
//   //vector<HMatrix<double>> SH = tt;
//   int of_t = t.get_offset(); int of_s= s.get_offset();
//   int sz_t = s.get_size(); int sz_s = s.get_size();
//   int n = SH.size()/2;
//   for(int k = 0 ; k< n; ++k){
//     HMatrix<double> bk = SH[k];
//     cout<< bk->get_target_cluster().get_offset()<< endl;}
//     int of_ = t.get_offset(); int of_s= s.get_offset();
//     int sz_t = s.get_size(); int sz_s = s.get_size();
  
//   return S;
//   }
// pair<vector<pair<Matrix<double>,int>>,pair< VirtualCluster,VirtualCluster>> Restrict (pair<vector<pair<Matrix<double>,int>>,pair< VirtualCluster,VirtualClsuter>> S, VirtualCluter& t , VirtualCluster& s){
//   vector<pair<Matrix<double>,int>> SR0 = S.first;
//   return (S,(t,s));
//}
// pair<pair<vector<Matrix<double>>,vector<int>>,vector<HMatrix<double>>> Restrict(pair<pair<vector<Matrix<double>>,vector<int>>,vector<HMatrix<double>>> S, const VirtualCluster &t, const VirtualCluster &s){
//   pair<vector<Matrix<double>>,vector<int>> SR0 , SR;
//   //vector<HMatrix<double>> SH = S.second;
//   SR0 = S.first;
//   SR = restrict_lr(SR0,t,s);
//   //auto tt = S.second;
//   //cout<< typeid(tt).name()<< endl;
//   //vector<HMatrix<double>> SH = tt;
//   int of_t = t.get_offset(); int of_s= s.get_offset();
//   int sz_t = s.get_size(); int sz_s = s.get_size();
//   // int n = SH.size()/2;
//   // for(int k = 0 ; k< n; ++k){
//   //   HMatrix<double> bk = SH[k];
//   //   cout<< bk->get_target_cluster().get_offset()<< endl;}
//     //int of_ = t.get_offset(); int of_s= s.get_offset();
//     //int sz_t = s.get_size(); int sz_s = s.get_size();
  
//   return S;
//   }
 // //on initialise Sh0 a 0
  // std::vector<pair<Block<double>,Block<double>>> Sh = S.second;
  // std::vector<pair<Block<double>,Block<double>>> Sh0;
  // //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
  // for (auto it= Sh.begin(); it!= Sh.end(); ++it){
  //   pair<Block<double>,Block<double>> HK = *it; Block<double> H=HK.first; Block<double> K = HK.second;
  //   Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
  //   int repere1 = 0;
  //   //on parcours les enfants de rho=r
  //   while(!(r.get_son_ptr(repere1) == nullptr)){
  //     //on cherche les restrictions de H et K parmi leur fils comme pour Sr
  //     int repereH =0;int repereK=0;Block<double> H0;Block<double> K0;
  //     int offset_r0 = r.get_son_ptr(repere).get_offset();
  //     int size_r0= r.get_son_ptr(repere).get_size();
  //     //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
  //     while(!(&H.get_son(repereH) == nullptr)){
  // 	Block<double> &Htemp = H.get_son(repereH);
  // 	int offset_t0 = Htemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Htemp.get_source_cluster().get_offset();
  // 	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
  // 	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
  // 	  H0 = Htemp;}
  // 	repereH +=1;}
  //     // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
  //     while(!(&K.get_son(repereK) == nullptr)){
  // 	Block<double> &Ktemp = K.get_son(repereK);
  // 	int offset_t0 = Ktemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Ktemp.get_source_cluster().get_offset();
  // 	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
  // 	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
  // 	  K0 = Ktemp;}
  // 	repereK +=1;}
  //     //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
  //     // Si c'est le cas le produit est forcément low rank
  // 	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
  // 	//On veut faire ABt = H0*K0------>push_back->Sr0
  // 	// trois possibilité: lr*lr, lr*hmat, hmat*lr
  // 	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
  // 	// Dans les trois cas on a un lr en résultat
  // 	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
  // 	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
  // 	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
  // 	 LowRankMatrix<double>& hr0 = H0.get_low_rank_block_data();
  // 	 LowRankMatrix<double>& kr0 = K0.get_low_rank_block_data();
  // 	 DenseBlockData<double>& hf0 = H0.get_dense_block_data();
  // 	 DenseBlockData<double>& kf0 = K0.get_dense_block_data();
  // 	// //cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
  // 	// if((hf0 == nullptr) and (kf0 == nullptr)){
  // 	//   LowRankMatrix<T> hk0 = h0*k0;}
  // 	// //cas 2 : hmat*lowrank
  // 	// else if( hr0 == nullptr ){
  // 	//   Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
  // 	//   int rank = U.nb_cols();
  // 	//   LowRankMatrix<T> hk0;
  // 	//   Matrix<T>* hk0U = hk0.Get_U();
  // 	//   Matrix<T>* hk0V = res.Get_V();
  // 	//   *hk0V = V;
  // 	//   Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
  // 	//   for ( int k =0; k< rank;++k){
  // 	//     vector<T> Uk = U.get_col(k);
  // 	//     vector<T> temp = H0*Uk;
  // 	//     for ( int l =0; l< U.nb_rows();++l){
  // 	//       mattemp(l,k)= Uk[l];
  // 	//     }
  // 	//   }
  // 	//   *hk0U = mattemp;
  // 	// }
  // 	// //cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
  // 	// // Sinon on code un produit vecteur mat xA = At x
  // 	// else if( kr0 ==nullptr){
	  
  // 	// }
	  
  //     }
  //     else{
  // 	pair<Block<double>,Block<double>> HK0; HK0.first = H0;HK0.second = K0;
  // 	Sh0.push_back(HK0);
  //     }
  //   }

//////////////////////////////////
// vector<Block<double>*> restrict( vector<Block<double>*> SR, const VirtualCluster& t , const VirtualCluster& s){
//   // // On fait la restriction SR
//   vector<Block<double>*> Sr0;
//   int oft_ref = t.get_offset() ; int ofs_ref = s.get_offset();
//   int szt_ref = t.get_size(); int szs_ref = s.get_size();
//   // copy(SR.begin(), SR.end(), back_inserter(Sr0));
//   int taille = SR.size();
//   int nr = 0;
//   for (int k =0 ;k < taille ; ++k){
//     Block<double>* bk = SR[k];
//      if( bk->IsAdmissible()){
//        int of_t = bk->get_target_cluster().get_offset() ; int sz_t = bk->get_target_cluster().get_size();
//        int of_s = bk->get_source_cluster().get_offset() ; int sz_s = bk->get_source_cluster().get_size();
       
//       if ( ( (of_t <= oft_ref) and (sz_t >= szs_ref)) and ((of_s <= ofs_ref) and (sz_s >= szs_ref) ) ){
//   	const LowRankMatrix<double>* lr =bk->get_low_rank_block_data();
//   	if(!(lr== nullptr)){
	  
//   	  Matrix<double> U = lr->Get_U();Matrix<double> V = lr->Get_V();
//     	 // const Matrix<double> U,Ures,V,Vres;
//     	 // U = lrk->Get_U(); V = lrk->Get_V();
//     	 // On veut faire la restriction de U selon les lignes et V les collonnes
//     	  Matrix <double> Uref ( szt_ref, U.nb_cols());
//     	  Matrix <double> Vref(V.nb_rows(), szs_ref);
// 	  cout<<"=============="<< endl;
// 	  cout<< V.nb_cols()<<','<<of_t -oft_ref<<','<< of_t<<','<< oft_ref<< endl;
// 	  cout<< sz_t<<','<<szs_ref<< endl;
// 	  for( int i =0; i< szt_ref; ++i){
// 	    for( int j =0; j< U.nb_cols(); ++j){
// 	      Uref(i,j) = U(i+of_t-oft_ref,j);}
//      }
// 	  for( int i =0; i< V.nb_rows(); ++i){
//       for (int j=0; j< szs_ref; ++j){
// 	//cout<< szs_ref<<","<<j+of_s-ofs_ref<<","<<V.nb_cols()<< endl;
// 	Vref(i,j) = V(i,j+of_s-ofs_ref);
//      }
//     }
//     // cout<<"============"<<endl;
//     // cout<< bk->nb_sons()<<endl ;
//     // bk->build_son(t,s);
//     // cout<< bk->nb_sons()<<endl;
//     // temp->get_low_rank_block_data()->Get_U() = Uref  ;
//     // temp->get_low_rank_block_data()->Get_V() = Vref;
//     // cout<< "$"<< endl;
//     // auto temp = bk->get_low_rank_block_data();
//     // // temp->set_U(Uref);
//     //   temp->Get_V() ;
//     //   Matrix<double>  vv =  temp->Get_V();
//     //   Matrix<double>* vvv = &vv;
//     //   //cout<< vv<< endl;
//     //   cout<<"================="<< endl;
//     //   cout<< vv.nb_rows()<<','<< vv.nb_cols()<< endl;
//     //   cout<< Vref.nb_rows()<<','<< Vref.nb_rows()<< endl;
//       //&vv = Vref;
//      //for (int k = 0; k < temp->Get_U().nb_rows(),++k){
//      // for( int l = 0 ; l< temp->Get_V.nbcols(),++k){
//      // if
//      // for( int k =0; k < temp->nb_rows();++k){
//      //    for ( int l =0; l< temp->nb_cols(); ++k){
//      //  	 double* val = &Uref(k,l);
//      //  	 temp->assign_U(k,l,*val);}
//      //  }
     
     
//      //UU=Uref ; VV=Vref;
//     // cout<< UU.nb_rows() <<endl;
    
//     //cout<< aa->nb_rows()<< endl;
//     //if (!(aa != nullptr)){
//       //  aa->Get_U();}
//       //  cout<< bk->nb_rows()<< endl;}
//    //  bk->Get_U() = Uref; bk->Get_V() = Vref;
//    //   bk->get_target_cluster() = t; bk->get_source_cluster()=s;
//    //   Block<double>* Bref ;
//    //   LowRankMatrix<double> Btemp = *Bref->get_low_rank_block_data();
//    //   Btemp.Get_U() = Uref;
//    //   Btemp.Get_V() = Vref;
//    //       Sr0.push_back(bk);
//    //   Matrix<double>* tempU = Bref.Get_U(); Matrix<double>* tempV= Bref.Get_V();
//    //   *tempU= Uref; *tempV= Vref
//    //   bk->build_son(t,s);
//    //  Block<double>& btemp = bk->get_son(0);
//     //Sr0.push_back(bk);
// 	}}
//      }
//    }

  

//    return SR;}
  
  
						    
////////////////////////////////////////////////////////////////////////////////////


//pair< vector<Block<double>>,vector<pair<Block<double>>,Block<double>>> restrict( pair< vector<Block<double>>,vector<pair<Block<double>>,Block<double>>> S, VirtualCluster &t, VirtualCluster &s){
  //On commence avec Sr
  // en faite on fera less troncatures a la fin et pour l'instant on fiat qu'en rajouter
  //   std::vector<Block<double>> Sr0;    int offset_t = t.get_offset(); int offset_s = s.get_offset(); int size_t = t.get_size(); int size_s = s.get_size();
  // for (auto it= Sr.begin(); it!= Sr.end(); ++it){
  //   // Block<double> &lr = *it;
  //   int repere = 0;
  //   // // ____________________//
  //   // //Pas sure que cette méthode marche mais sinon il faut mettre la matrice en full, faire sa restriction et refaire ACA;
  //   // // on peut aussi juste faire une boucle for k < lr.get_nb_son()
  //   // //____________________//
  //   // // On récupère la réstriction a tau sigma en vérifiant sur chaques blocs des get_son() si leurs clusters sont les mêmes
  //   // while(!(&lr.get_son(repere) ==nullptr)){
  //   //   Block<double> &lr_son = lr.get_son(repere);
  //   //   int offset_t0=lr_son.get_target_cluster().get_offset(); int offset_s0= lr_son.get_source_cluster().get_offset();
  //   //   int size_t0= lr_son.get_target_cluster().get_size(); int size_s0 = lr_son.get_source_cluster().get_size();
  //   //   // on vérifie que les cluster sont les mêmes <=> même sze et offset
  //   //   if( ((offset_t == offset_t0) and (offset_s == offset_s0))  and ((size_t0 == size_t) and (size_s0 == size_s))){
  //   // 	Sr0.push_back(lr_son);}
  //     repere +=1;
  //   }
  // vector<Block<double>> Sr0 = S.first;
  
  // //on initialise Sh0 a 0
  // std::vector<pair<Block<double>,Block<double>>> Sh = S.second;
  // std::vector<pair<Block<double>,Block<double>>> Sh0;
  // //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
  // for (auto it= Sh.begin(); it!= Sh.end(); ++it){
  //   pair<Block<double>,Block<double>> HK = *it; Block<double> H=HK.first; Block<double> K = HK.second;
  //   Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
  //   int repere1 = 0;
  //   //on parcours les enfants de rho=r
  //   while(!(r.get_son_ptr(repere1) == nullptr)){
  //     //on cherche les restrictions de H et K parmi leur fils comme pour Sr
  //     int repereH =0;int repereK=0;Block<double> H0;Block<double> K0;
  //     int offset_r0 = r.get_son_ptr(repere).get_offset();
  //     int size_r0= r.get_son_ptr(repere).get_size();
  //     //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
  //     while(!(&H.get_son(repereH) == nullptr)){
  // 	Block<double> &Htemp = H.get_son(repereH);
  // 	int offset_t0 = Htemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Htemp.get_source_cluster().get_offset();
  // 	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
  // 	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
  // 	  H0 = Htemp;}
  // 	repereH +=1;}
  //     // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
  //     while(!(&K.get_son(repereK) == nullptr)){
  // 	Block<double> &Ktemp = K.get_son(repereK);
  // 	int offset_t0 = Ktemp.get_target_cluster().get_offset();
  // 	int offset_s0 = Ktemp.get_source_cluster().get_offset();
  // 	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
  // 	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
  // 	  K0 = Ktemp;}
  // 	repereK +=1;}
  //     //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
  //     // Si c'est le cas le produit est forcément low rank
  // 	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
  // 	//On veut faire ABt = H0*K0------>push_back->Sr0
  // 	// trois possibilité: lr*lr, lr*hmat, hmat*lr
  // 	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
  // 	// Dans les trois cas on a un lr en résultat
  // 	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
  // 	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
  // 	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
  // 	 LowRankMatrix<double>& hr0 = H0.get_low_rank_block_data();
  // 	 LowRankMatrix<double>& kr0 = K0.get_low_rank_block_data();
  // 	 DenseBlockData<double>& hf0 = H0.get_dense_block_data();
  // 	 DenseBlockData<double>& kf0 = K0.get_dense_block_data();
  // 	// //cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
  // 	// if((hf0 == nullptr) and (kf0 == nullptr)){
  // 	//   LowRankMatrix<T> hk0 = h0*k0;}
  // 	// //cas 2 : hmat*lowrank
  // 	// else if( hr0 == nullptr ){
  // 	//   Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
  // 	//   int rank = U.nb_cols();
  // 	//   LowRankMatrix<T> hk0;
  // 	//   Matrix<T>* hk0U = hk0.Get_U();
  // 	//   Matrix<T>* hk0V = res.Get_V();
  // 	//   *hk0V = V;
  // 	//   Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
  // 	//   for ( int k =0; k< rank;++k){
  // 	//     vector<T> Uk = U.get_col(k);
  // 	//     vector<T> temp = H0*Uk;
  // 	//     for ( int l =0; l< U.nb_rows();++l){
  // 	//       mattemp(l,k)= Uk[l];
  // 	//     }
  // 	//   }
  // 	//   *hk0U = mattemp;
  // 	// }
  // 	// //cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
  // 	// // Sinon on code un produit vecteur mat xA = At x
  // 	// else if( kr0 ==nullptr){
	  
  // 	// }
	  
  //     }
  //     else{
  // 	pair<Block<double>,Block<double>> HK0; HK0.first = H0;HK0.second = K0;
  // 	Sh0.push_back(HK0);
  //     }
  //   }
  // return Sr0;
  // }
  //  }

























    
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
    double epsilon = 0.00001;
    double eta     = 0.5;

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
    //Aller on test
    vector<Block<double>*> Bl = HA.get_ComputedBlock();
    //vector<pair<int,Matrix<double>>,pair<int,Matrix<double>>> SR;
    //pair<vector<pair<int,Matrix<double>>>,vector<pair<int,Matrix<double>>>> SR;
    pair<vector<Matrix<double>>,vector<int>> SR;
    vector<Matrix<double>> sr; vector<int> of;
    for (int k = 0; k < Bl.size(); ++k){
      Block<double>* bk = Bl[k];
      int ofu = bk->get_target_cluster().get_offset(); int szu = bk->get_target_cluster().get_size();
      int ofv = bk->get_target_cluster().get_offset(); int szv = bk->get_source_cluster().get_size();
      //pair<int,int>
      pair<int, Matrix<double>> Uref,Vref;
      if( !( bk->get_low_rank_block_data()== nullptr) ){
	Matrix<double> U = bk->get_low_rank_block_data()->Get_U();
	Matrix<double> V = bk->get_low_rank_block_data()->Get_V();
	sr.push_back(U);
	sr.push_back(V);
	of.push_back(ofu); of.push_back(ofv);
      //Uref(ofu,U); Vref(ofv,V);
      //pair<pair<int,Matrix<double>>,pair<int,Matrix<double>>> temp (Uref,Vref);
      //SR.push_back(temp);
       
    }
    }
    SR.first = sr; SR.second = of;
    cout<< sr.size()<<"'"<<of.size()<< endl;
    pair<vector<Matrix<double>>,vector<int>> Restr =restrict_lr (SR,*t,*t);
    vector<Matrix<double>> aa = Restr.first;
    // for(int ii = 0 ; ii< aa.size(); ++ii){
    //   Matrix<double> M = aa[ii];
    //   Matrix<double> MM = sr[ii];
    //   cout<< M.nb_rows()<<','<<M.nb_cols()<<endl;
    //   cout<< MM.nb_rows()<<','<< MM.nb_cols()<< endl;
    //   cout<<"===================="<<endl;}
    HMatrix<double>* HH =&HA;
     HH->get_BlockTree() ;
     cout<< HH->get_target_cluster()->get_size()<<endl;


     //vector<const VirtualCluster&> pp;pp.push_back(t);
     //cout <<pp.size()<<endl;
     cout<< t->get_nb_sons()<<endl;
     cout<< t->get_son(0).get_son(0).get_nb_sons()<<endl;
     VirtualCluster& i= t->get_son(0);
     cout<<"$$$$$$$$$$$$$$$"<< endl;
     cout<< i.get_size()<<endl;
     vector<VirtualCluster*> vec ;
     vec.push_back(&i);
     cout<< vec.size()<<endl;
     cout<< vec[0]->get_size()<< endl;
     vec.push_back(&i);
     //pair< pair<vector<Matrix<double>>,vector<int>> , vector<VirtualCluster>&> lll;
     //lll.first = Restr;
     //lll.second = vec;
     pair< pair<vector<Matrix<double>>,vector<int>> , vector<pair<pair<int,int>,pair<int,int>>>> lll;
     lll.first = Restr;

     Block<double>* bl0 = Bl[2];
     cout<< bl0->get_target_cluster().get_size()<< endl;
     cout<< bl0->get_root()->get_target_cluster().get_size()<< endl;
     cout<< bl0->get_root()->get_son(0).get_son(0).get_target_cluster().get_size()<< endl;
     cout<< bl0->get_root()->get_son(0).get_son(0).get_son(0).get_target_cluster().get_size()<< endl;
     Block<double>* rtt = bl0->get_root();
     //while(!(rtt.get_son(0) ==nullptr)){
     //cout<<rtt.get_target_cluster().get_size();<< endl;
     //rtt= rtt.get_son(0);}
     int ttn = rtt->nb_sons();
     cout<< ttn << endl;

     //
     //Test Restrict//
     ////////////////
     cout<< "========================"<< endl;
     cout<< " Test Restrict"<< endl;
     cout<<"========================"<<endl;
     Block<double>* Rt = bl0->get_root();vector<Block<double>*> SH; SH.push_back(Rt); SH.push_back(Rt);
     pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> T0; T0.first= Restr;T0.second = SH;
     pair<pair < vector<Matrix<double>> ,vector<int>>,vector<Block<double>*>> TEST = Restrict(Restr,SH,*t,*t);
     int rr = 0;
     cout<< Rt->nb_sons()<<endl;
     cout<< Rt->get_son(0).nb_sons()<<endl;
     //Block<double>* tempp = new Block<double>;
     int nb = Rt->get_son(0).nb_sons();
     cout<< nb<< endl;
     Block<double>& tempp = Rt->get_son(0);
     cout<< tempp.nb_sons()<< endl;
     cout<< bl0->nb_sons()<< endl;
     Block<double>* bk = &tempp;
     int bb = 0;
     while( bb<3){
       
       cout<< bk->get_target_cluster().get_size()<<endl;
       Block<double>& bbk = bk->get_son(0);
       bk =&bbk;
       
       //Block<double>bk =tp;
       //tempp= tp;
       bb+=1;}
       
     cout<< "test prod"<<endl;
     cout<<"===================="<< endl;
     vector<double> xx (Rt->get_source_cluster().get_size(),1);
     vector<double> xxx(Rt->get_source_cluster().get_size(),1);
     vector<double> yy (Rt->get_target_cluster().get_size(),0);
     Block<double>* ttemp = &Rt->get_son(0);
     vector<double> yp = Prod(ttemp,xx,yy);
     vector<double> yref =HA*xx;
     vector<double> ytest = A*xxx;
     double tt=0;
     cout<< yp[10]<<','<<ytest[10]<<endl;
     int ss =0;
     //for (int ut =0 ; ut< ytest.size();++ut){
     //ss+=(ytest[ut]*ytest[ut]);}

          for (int ut =0 ; ut< ytest.size();++ut){
        for(int us = 0; us< A.nb_cols();++us){
        	 ss+=(A.get_coef(ut,us)*A.get_coef(ut,us));}}

     for (int k =0; k < ytest.size(); ++k){
       tt+= ((yp[k]-ytest[k])*(yp[k]-ytest[k]));}
     cout<< sqrt(tt/ss)<< endl;
     cout<<NORM(yp,yref)/sqrt(ss)<<endl;
     cout<<sqrt(100)<<endl;

     // int  zer = 0;
     vector<double> yz(Rt->get_target_cluster().get_size(),0);
     vector<double>* yyz= &yz;
     // vector<int> repp;
     // Produit(ttemp,yyz,xxx);
     //vector<double> val = *yyz;
     // cout<< NORM(val,ytest)/sqrt(ss)<<endl;
     // for(int k =0; k< 20; ++k){
     //   cout<< val[k]<<','<<ytest[k]<<endl;}
     // cout<<NORM(yref,ytest)/sqrt(ss)<<endl;
     // cout<< Rt->get_target_cluster().get_offset()<<endl;
     // cout<< Rt->get_target_cluster().get_size()<<endl;
     // cout<< Rt->get_source_cluster().get_offset()<<endl;
     // cout<< Rt->get_source_cluster().get_size()<<endl;
     cout<<";;;;;;;;;;;;;;;;;;;;;;;;;"<<endl;
     cout<< Rt->get_tasks().size()<<endl;
     cout<< Rt->get_local_tasks().size()<<endl;
     vector<Block<double>*> St = Rt->get_tasks();
     int nn=0;
     int nt = 0; int nr = 0;
      // for (int k =0 ; k< St.size(); ++k){
     //    Block<double>* s0 = St[k];
     // 	vector<Block<double>*> sss = s0->get_tasks();
     // 	if((s0->get_dense_block_data()==nullptr) and (s0->get_block_data()==nullptr)){nn+=1;}
     // //   cout<< s0->get_target_cluster().get_offset()<<","<<s0->get_target_cluster().get_size()<<endl;
     // //   cout<< s0->get_source_cluster().get_offset()<<","<<s0->get_source_cluster().get_size()<< endl;}
     //    if(s0->get_block_data()==nullptr){nt +=1;}
     //    if(s0->get_dense_block_data()==nullptr){nr +=1;}
     // //   cout<< s0->get_target_cluster().get_offset()<<","<<s0->get_target_cluster().get_size()<<endl;
     // //   cout<< s0->get_source_cluster().get_offset()<<","<<s0->get_source_cluster().get_size()<< endl;}
     // 	cout<< sss.size()<<endl;
     
     // }
      cout<<nn<<","<<nt<<","<<nr<< endl;
      cout<< Rt->get_local_tasks().size()<< endl;//Block<double>* s = St[0];cout<< s->get_local_diagonal_block().get_target_cluster().get_size()<<endl;
     // cout<< Rt->get_local_diagonal_block().get_target_cluster().get_size()<< endl;
      vector<Block<double>*> lst = HA.get_ComputedBlock();
      // cout<<lst.size()<< endl;
      // double nnx= 0; double nny =0;
      // double nnr = 0;
      // int  nlr=0; double testn=0;double test1=0;
      // double testt = 0; double testr = 0;
      // int cpt =0;
      // for(auto it = lst.begin(); it != lst.end();++it){
      // 	Block<double>* l = *it;
      // 	int of_t = l->get_target_cluster().get_offset(); int sz_t= l->get_target_cluster().get_size();
      // 	int of_s = l->get_source_cluster().get_offset(); int sz_s = l->get_source_cluster().get_size();
      // 	vector<double> xtest(sz_s,1); vector<double> ytest(sz_t,0);
      // 	if(l->get_low_rank_block_data()==nullptr){
      // 	  cout<<cpt<<",";
      // 	Matrix<double> mm(sz_t,sz_s);
      // 	//if(l->get_block_data()==nullptr){
      // 	  //cout<<"!"<<r12<<endl;}
      // 	//r12+=1;
      // 	 for (int i = 0; i< sz_t; ++i){
      // 	     for(int j = 0; j< sz_s; ++j){
      // 	       mm(i,j)=A.get_coef(i+of_t,j+of_s);}}
      // 	 l->get_block_data()->add_mvprod_row_major(xtest.data(),ytest.data(),1,'T','N');
      // 	 vector<double> xtest2(sz_t,1);
      // 	 vector<double> yres = mm*xtest2;
      // 	 //cout<<NORM(ytest,yres)<< endl;
      // 	 for (int k=0; k < sz_t;++k){
      // 	   nnx+= ytest[k]*ytest[k]; nny+= yres[k]*yres[k];}
      // 	 nnr+= NORM(ytest,yres);}
      // 	//if(l->get_block_data()==nullptr){cout<<"merde"<<endl;}
      // 	else{
      // 	  Matrix<double> U= l->get_low_rank_block_data()->Get_U(); Matrix<double> V = l->get_low_rank_block_data()->Get_V();
      // 	  vector<double> xy (sz_s,1);
      // 	  vector<double> yx  (sz_t,0);
      // 	  l->get_low_rank_block_data()->add_mvprod_row_major(xy.data(),yx.data(),1,'T','N');
      // 	  vector<double> yyx = U*(V*xy);
      // 	  Matrix<double> mm(sz_t,sz_s);
      // 	//if(l->get_block_data()==nullptr){
      // 	  //cout<<"!"<<r12<<endl;}
      // 	//r12+=1;
      // 	 for (int i = 0; i< sz_t; ++i){
      // 	     for(int j = 0; j< sz_s; ++j){
      // 	       mm(i,j)=A.get_coef(i+of_t,j+of_s);}}
	
      // 	 vector<double> xtest2(sz_t,1);
      // 	 vector<double> yres = mm*xtest2;
      // 	 //cout<<NORM(ytest,yres)<< endl;
      // 	 for (int k=0; k < sz_t;++k){
      // 	   testt+= yx[k]*yx[k]; testr+= yres[k]*yres[k];}
      // 	 testn+= NORM(yx,yres);test1+=NORM(yyx,yres);}
      // 	cpt+=1;
      // }
      // cout<< sqrt(nnr)<<","<< sqrt(nnx)<<","<<sqrt(nny)<< endl;
      // cout<<testn<<","<<test1<<','<<sqrt(testt)<<","<<sqrt(testr)<<endl;
      double sz1=0; double sz2=0;int cpt=0;
      for (int mom = 0; mom<lst.size(); ++ mom){
	if(lst[mom]->get_low_rank_block_data() == nullptr){
	  cpt+=1;
	//const LowRankMatrix<double>* lr0 = lst[18151]->get_low_rank_block_data();
	int sz_t =lst[mom]->get_target_cluster().get_size();
        int sz_s =lst[mom]->get_target_cluster().get_size();
	  Matrix<double> m1(sz_t,sz_s);
      for(int i=0; i< sz_t; ++i){
      	for(int j =0; j< sz_s;++j){
      	  m1(i,j)=A.get_coef(i+lst[mom]->get_target_cluster().get_offset(),j+lst[mom]->get_source_cluster().get_offset());
      	}
      }
      vector<double> x1(sz_s,1);
      vector<double> x2(sz_s,1);
      vector<double> y1(sz_t,0);
      vector<double> y2(sz_t,0);
      //cout<<"++++++++"<<endl;
      lst[mom]->get_block_data()->add_mvprod_row_major(x1.data(),y1.data(),1,'N','N');
      //cout<< NORM(m1*x1,y1)/sqrt(lst[mom]->get_target_cluster().get_size())<<endl;
      lst[mom]->get_block_data()->add_mvprod_row_major(x2.data(),y2.data(),1,'T','N');
      sz1+= NORM(y1,m1*x1)*NORM(y1,m1*x1); sz2+= NORM(y2,m1*x2)*NORM(y2,m1*x2);
      //cout<< NORM(m1*x1,y2)/sqrt(lst[mom]->get_target_cluster().get_size())<<endl;
      //cout<< lst[mom]->get_target_cluster().get_offset()<<","<<lst[mom]->get_target_cluster().get_size()<<endl;
      // Matrix<double> U = lr0->Get_U(); Matrix<double> V =lr0->Get_V();Matrix<double> W = U*V;
      // for(int k =0; k< U.nb_rows();++k){
      // 	for(int l =0; l< V.nb_cols();++l){
      // 	  cout<< m1(k,l)<<"!"<<W(k,l)<< '\t';}
      // 	cout<<endl;}
	}}
      cout<<sqrt(sz1/cpt)<<','<<sqrt(sz2/cpt)<<endl;




















     
     //vector<double> yyref = Prod(ttemp,xx,yy);
     // for (int k =0; k < ytest.size(); ++k){
     //   tt+= (ytest[k]-yyref[k])*(ytest[k]-yyref[k]);}

     // cout<< sqrt(tt/ss)<< endl;

     //Block<double>* r1 = new Block<double>;
     //r1 =Rt->get_son(0);
     //auto it = Rt->get_son(0);
     //for ( int k = 0; k<5; ++k){
       
     //Block<double> rt = Rt->get_son(0);
     //Block<double>* Rtt;
     //Rtt = Rt;
     //Block<double> tempp = Rt->get_son(0);
     //     while (rr<3){
       //       cout<< tempp.get_target_cluster().gt_size()<<enl;
       //	  tempp = tempp.get_son(0);
     //   Block<doub
     //   Block<double> ttt = rtt->get_son(0);
     //   cout<<"+++++++"<< endl;
     //   cout<<Rt->nb_sons()<<endl;
     // //cout<<Rt->get_target_cluster().get_size() <<endl;
	  //	rr+=1;
     //   rtt = &ttt;
     //   //cout<< rtt->get_target_cluster().get_size()<< endl;
	//}
     
     //pair<int,int
     //Restrict( (Restr,vec),*t,*t);
     //cout<< vec.size()<< endl;
     //vec.push_back(i);
     //cout<< vec.size()<< endl;
     //vector< VirtualCluster> vec;
     //vec.push_back(t->get_son(0).get_son(0));
     //Cluster t1= t->get_son(1);
     //VirtualCluster &t1 = t.get_son(1);
     //cout<< t1.get_size()<< endl;

     //Block<double> B1 = HH->get_BlockTree()->get_son(1);
     //cout<< B1.get_target_cluster().get_size()<< endl;
    
    
    //HA.get_BlockTree();

    // pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> Test;
    // pair<Block<double>*,Block<double>*> temp;
    // temp.first = HA; temp.second = HA;
    //Test.fisrt = Restrict(// Block<double>* b1 = Bl[10];
    // cout<< b1->get_source_cluster().get_offset()<< endl;
    // VirtualCluster vv = b1->get_source_cluster();
    // cout< "+++++++++++++++++"<< endl;
    // for (int k =0; k < vv.size(); ++k){
    //   cout << vv[k<<endl;}
    
    // HA.save_plot(outputpath + "/smallest_example_plot");
    // HA.get_target_cluster()->save_geometry(p.data(), outputpath + "/smallest_example_cluster", {1, 2, 3});
    // std::cout << outputpath + "/smallest_example_plot" << std::endl;
    // std::cout << Frobenius_absolute_error(HA, A) / A.norm() << std::endl;
    // std::cout << norm2(A * x - result) / norm2(A * x) << std::endl;


    // On rajoute nos test.


    // On veut modifier les blocs de L =HK

    // On copie L ( j'arrive pas a l'initialiser a 0
    // HMatrix<double> Res(t, t, epsilon, eta, 'S', 'U');
    // Res.build(A, p.data());
    // vector<Block<double>*> Data = HA.get_ComputedBlock();
    // vector<Block<double>*> DataRes = Res.get_ComputedBlock();
    // // On va classer les low rank et pas low rank


//     vector<Block<complex<double>>> ll,fl,llR,flR;

//     for(auto it = Data.begin(); it!= Data.end(); ++it){
//       int ofs_t =  it->get_target_cluster.get_offset(); int ofs_s  = it->get_source_cluster.get_offset();
//       int sz_t =  it->get_target_cluster.get_size(); int sz_s  = it->get_source_size.get();
//       if( 
// fonction restrict

// test pour la restriction low rank

// d'abord on fait la liste des low rank

    // vector<Block<double>*> liste = HA.get_ComputedBlock();
    // vector<Block<double>*> ll;
    // int taille = liste.size();
    // for (int k =0 ; k < taille; ++k){
    //   if (liste[k]->IsAdmissible()){
    // 	ll.push_back(liste[k]);
    // 	  Block<double>* kk = liste[k];
    // 	  const LowRankMatrix<double>* lrr = kk->get_low_rank_block_data();
    // 	  if (!(lrr==nullptr)){
    // 	    lrr->Get_U();
    //   }
    // }
    // }
    // cout<<"taille"<<ll.size()<<endl;
    // for(int k =0 ; k< 10; ++k){
    //   cout<< ll[k]->IsAdmissible()<< endl;}
    // //vector<Block<double>*> test = restrict(ll,*t,*t);
    // cout<<"================"<< endl;
    // // t->get_global_perm();
    // // for(auto it = t->get_global_perm().begin(); it!= t->get_global_perm().end(); ++it){
    // //   cout<< *it<< endl;}
    // // Block<double>* bltest = test[10];
    // // const LowRankMatrix<double>* lrtest = bltest->get_low_rank_block_data();
    // // Block<double>* bltest2 = liste[10];
    // // const LowRankMtarix<double> lrtest2 = bltest2->get_lowrank_bock_data();
    // // cout<< lr


    // Block<double>* tt = ll[20];
    // //Block<double> tep = &tt;
    //  // auto tete = tt->get_low_rank_block_data();
    //  //cout<< tete->get_U(1,1) <<endl;
    // // // cout<<"t et s" << endl;
    // // // cout<< tt->get_target_cluster().get_size()<< ','<< tt->get_source_cluster().get_size()<< endl;
    // // // cout<< tt->get_target_cluster().get_offset()<< ','<< tt->get_source_cluster().get_offset();
    //  Matrix<double>* M1,M2,M3;
    //  cout<< tt->IsAdmissible()<< endl;
    // const LowRankMatrix<double>* lr =tt->get_low_rank_block_data();
    // 	if(!(lr== nullptr)){

    // 	  //Matrix<double> UU =&lr->Get_U();
    // 	  Matrix<double> U = lr->Get_U();Matrix<double> V = lr->Get_V();
    // 	  cout<<U(0,0)<< endl;
    // 	  Matrix<double>* UU =&U;
    // 	  double a = U(0,0);
    // 	  Matrix<double> AA (10,10);
    // 	  for(int k =0; k <U.nb_rows() ; ++k){
    // 	    for (int l =0; l<U.nb_cols(); ++l){
    // 	    AA(k,l)=99;
    // 	  }
    // 	}
    // 	  cout<<U(0,0)<< endl;
    // 	}
	  //cout<< lr->Get_U()(0,0)<< endl;
	  //Block<double>& temp = *tt;
	  //temp.get_low_rank_block_data()->Get_U()(0,0)= 9;
	  //double* a = 9.8;
	  //temp.get_low_rank_block_data()->assign_U(0,0,a);
	  //cout<< temp.get_low_rank_block_data()->Get_U()(0,0)<< endl;
	  //LowRankMatrix<double> testM( AA, V);
	  
	  
	  // LowRankMatrix<double> ( tt, A,lr,  10);
	//   Matrix<double> MM1 = tt->get_low_rank_block_data()->Get_U();
	//   cout<< MM1(0,0)<< endl;
    // 	      vector<Block<double>*> liste1 = HA.get_ComputedBlock();
    // vector<Block<double>*> ll1;
    // for (int k =0 ; k < taille; ++k){
    //   if (liste1[k]->IsAdmissible()){
    // 	ll1.push_back(liste[k]);
    // 	  Block<double>* kk = liste[k];
    // 	  const LowRankMatrix<double>* lrr = kk->get_low_rank_block_data();
    // 	  if (!(lrr==nullptr)){
    // 	    lrr->Get_U();
    //   }
    // }
    // }
    // Block<double>* UUU = ll1[20];
    //   const LowRankMatrix<double>* lr1 =tt->get_low_rank_block_data();
    // 	if(!(lr== nullptr)){
	  
    // 	  Matrix<double> U0 = lr1->Get_U();Matrix<double> V = lr->Get_V();
    // 	  cout<<U0(0,0)<< endl;
    // 	  //Matrix<double>* UUU =&U;
    // 	  Matrix<double> AA (10,10);
    // 	  for(int k =0; k <10 ; ++k){
    // 	    for (int l =0; l<10; ++l){
    // 	    AA(k,l)=99;
    // 	  }
    // 	}
    // 	}
	  //*UUU = AA;
	  // Matrix<double>* U0 = &U;
	  // U(0,0) =5;
	  // cout<< U(0,0) << endl;
	  // Matrix<double> MM = *U0;
	  // double *a;
	  // double  b = 5.8;
	  // MM(0,0) = 5.8;
	  // a=&b;
	  // cout<< U(0,0)<< endl;
	  

	  //LowRankMatrix<double> lrm( tt,A,lr,lr->get_xr(),lr->get_xc(),lr->rank_of());
	  
	  // for(int k = 0 ; k < U.nb_rows(); ++k){
	  //   for(int l =0 ; l< U.nb_cols(); ++l){
	  //     cout<< U(k,l)<<'\t';}
	  //   cout<< endl;}
	  // cout<< U(0,0)<< endl;
	  // double *u1;
	  // double u = U(0,0);double a = 5.9;
	  // u1=&u;
	  // u1=5.8;
	  // cout<< U(0,0)<< endl;
	  //}

    // vector<Block<double>*> liste1 = HA.get_ComputedBlock();
    // vector<Block<double>*> ll1;
    // for (int k =0 ; k < taille; ++k){
    //   if (liste1[k]->IsAdmissible()){
    // 	ll1.push_back(liste1[k]);

    // }
    // }
    //  const LowRankMatrix<double>* lr1 =tt->get_low_rank_block_data();
    // 	if(!(lr1== nullptr)){
	  
    // 	  Matrix<double> U = lr1->Get_U();Matrix<double> V = lr1->Get_V();
    // 	  // for(int k = 0 ; k < U.nb_rows(); ++k){
    // 	  //   for(int l =0 ; l< U.nb_cols(); ++l){
    // 	  //     cout<< U(k,l)<<'\t';}
    // 	  //   cout<< endl;}
    // 	  cout<< U(0,0)<< endl;
    // 	}

    //int nb = tt->get_low_rank_block_data()->rank_of();
    //cout<< nb<< endl;
    //*M1 = lrr->Get_U();
    // // M2 = lrr.Get_V();
    // // M3 = M1*M2;
    // // for (int k =0; k<  tt->get_target_cluster().get_size(); ++k){
    // //   for ( int l =0 ;  l< tt->get_source_cluster().get_size();++l){
    // // 	cout<< A(k+ tt->get_target_cluster().get_offset(),l+ tt->get_source_cluster().get_offset())-M3(k,l)<< '\t';}
    // // 	      cout<< endl;}
    // for (int k =0; k < test.size(); ++k){
    //   Block<double>* tt = test[k];
    //    cout<< tt->get_target_cluster().get_size()<< ','<< tt->get_source_cluster().get_size()<< endl;
    //    cout<< tt->get_target_cluster().get_offset()<< ','<< tt->get_source_cluster().get_offset()<<endl;}

    // Block<double>* tt = ll[20];
    // cout<< tt->get_target_cluster().get_size()<< ','<< tt->get_source_cluster().get_size()<< endl;
    // cout<< tt->get_target_cluster().get_offset()<< ','<< tt->get_source_cluster().get_offset()<<endl;
    // cout<<  tt->get_target_cluster().get_nb_sons()<< endl;
    // cout<< tt->get_target_cluster().get_son(0).get_size()<< ','<< tt->get_target_cluster().get_son(0).get_size()<< endl;
    // cout<< tt->get_target_cluster().get_son(0).get_offset()<< ','<< tt->get_source_cluster().get_son(0).get_offset()<<endl;
    // tt->build_son(tt->get_target_cluster().get_son(0),tt->get_target_cluster().get_son(0));
    //  cout<<  tt->get_target_cluster().get_nb_sons()<< endl;
   
    // for( int k=0; k < taille ; ++k){
    //   Block<double>* ll1 = ll[k];
    //   Block<double>* ll2 = test[k];
    //   int n = 0;
    //   Matrix<double> res1 = ll1->get_low_rank_block_data()->Get_U() * ll1->get_low_rank_block_data()->Get_V();
    //   Matrix<double> res2 = ll2->get_low_rank_block_data()->Get_U() * ll2->get_low_rank_block_data()->Get_V();
    //   // cout<<"taill"<< ll1->get_target_cluster().get_size()<< endl;
    //   // for ( int k =0 ; k < ll1->get_target_cluster().get_size(); ++k){
    //   // 	for(int l =0; l< ll2->get_source_cluster().get_size(); ++l){
    //   // 	  cout<<  sqrt(res1(k,l)*res1(k,l)-res2(k,l)*res2(k,l))<< endl;
    //   // 	}
    //   // }
    //   // cout << n << endl;
      
    // }
      cout<"whaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaat"<<endl;
      MPI_Finalize();
}
