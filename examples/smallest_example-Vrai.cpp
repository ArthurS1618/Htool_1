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
};

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
    HMatrix<double> HA(t, t, epsilon, eta, 'S', 'U');
    HA.build(A, p.data());
    result = HA * x;

    // Output
    HA.print_infos();
    HA.save_plot(outputpath + "/smallest_example_plot");
    HA.get_target_cluster()->save_geometry(p.data(), outputpath + "/smallest_example_cluster", {1, 2, 3});
    std::cout << outputpath + "/smallest_example_plot" << std::endl;
    std::cout << Frobenius_absolute_error(HA, A) / A.norm() << std::endl;
    std::cout << norm2(A * x - result) / norm2(A * x) << std::endl;


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

template < typename T>
  T= double;

std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> restrict(std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> S, Cluster t, Cluster s){
  std::vector<Block<T>> Sr = S.first;
  //il faut faire la restriction des low rank a t,s
  std::vector<Block<T>> Sr0;    int offset_t = t.get_offset(); int offset_s = s.get_offset(); int size_t = t.get_size(); int size_s = s.get_size();
  for (auto it= Sr.begin(); it!= Sr.end(); ++it){
    Block<T> lr = *it;
    int repere = 0;
    // ____________________//
    //Pas sure que cette méthode marche mais sinon il faut mettre la matrice en full, faire sa restriction et refaire ACA;
    // on peut aussi juste faire une boucle for k < lr.get_nb_son()
    //____________________//
    // On récupère la réstriction a tau sigma en vérifiant sur chaques blocs des get_son() si leurs clusters sont les mêmes
    while(!(&lr.get_son(repere) ==nullptr)){
      Block<T> &lr_son = lr.get_son(repere);
      int offset_t0=lr_son.get_target_cluster().get_offset(); int offset_s0= lr_son.get_source_cluster.get_offset();
      int size_t0= lr_son.get_target_cluster().get_size(); int size_s0 = lr_son.get_source_cluster().get_size();
      // on vérifie que les cluster sont les mêmes <=> même sze et offset
      if( ((offset_t == offset_t0) and (offset_s == offset_s0))  and ((size_t0 == size_t) and (size_s0 == size_s))){
	Sr0.push_back(lr_son);}
      repere +=1;
    }
  }
  //on initialise Sh0 a 0
  std::vector<pair<Block<T>,Block<T>>> Sh = S.second;
  std::vector<pair<Block<T>,Block<T>>> Sh0;
  //On parcours les éléments de Sh, on fait la restriction et on avise selon les cas
  for (auto it= Sh.begin(); it!= Sh.end(); ++it){
    pair<Block<T>,Block<T>> HK = *it; Block<T> H=HK.first; Block<T> K = HK.second;
    Cluster r =  H.get_source_cluster();// aussi égale a K.get_target_cluster()
    int repere1 = 0;
    //on parcours les enfants de rho=r
    while(!(r.get_son_ptr(repere1) == nullptr)){
      //on cherche les restrictions de H et K parmi leur fils comme pour Sr
      int repereH =0;int repereK=0;Block<T> H0;Block<T> K0;
      int offset_r0 = r.get_son_ptr(repere).get_offset();
      int size_r0= r.get_son_ptr(repere).get_size();
      //on parcours les fils de h et on regarde si son cluster source  est bien le même que rho= r0 et son target = tau=t= target de la restriction S(t,s)
      while(!(&H.get_son(repereH) == nullptr)){
	Block<T> &Htemp = H.get_son(repereH);
	int offset_t0 = Htemp.get_target_cluster().get_offset();
	int offset_s0 = Htemp.get_source_cluster().get_offset();
	int size_t0=Htemp.get_target_cluster().get_size(); int s0= Htemp.get_cluster_source().get_size();
	if( (((offset_t0 == offset_t) and (size_t0 == size_t)) and ((offset_s0 == offset_r0) and (size_s0== size_s)){
	  H0 = Htemp;}
	repereH +=1;}
      // pareil avec K mais cette fois on vérifie que K.target= r et K.source = s = restriction S(t,s).source 
      while(!(&K.get_son(repereK) == nullptr)){
	Block<T> &Ktemp = K.get_son(repereK);
	int offset_t0 = Ktemp.get_target_cluster().get_offset();
	int offset_s0 = Ktemp.get_source_cluster().get_offset();
	int size_t0 Ktemp.get_source_cluster().get_size(); int size_s0 = Ktemp.get_source_cluster().size();
	if( ((offset_t0 == offset_r0) and (size_t0 == size_r0)) and ((offset_s0 == offset_s) and (size_s0 == size_s))){
	  K0 = Ktemp;}
	repereK +=1;}
      //On a les restrictions maintenant on regarde si un des deux est une feuille admissible
      // Si c'est le cas le produit est forcément low rank
	  if (&H0.IsAdmissible() or &K0.IsAdmissible()){
	//On veut faire ABt = H0*K0------>push_back->Sr0
	// trois possibilité: lr*lr, lr*hmat, hmat*lr
	//les deux dernier son équivalent il suffit d'apliquer un des deux a la transposé
	// Dans les trois cas on a un lr en résultat
	//TO DO: implémenter le prduit low_low_rank/Hmat-----> class mylrmat() pour les test sans toucher a lrmat ?
	// un seul des lr_data , dene_data est different de ptrnull si c'est une feuille (admissible ou pas) et les deux ptrnull sinon
	// ainsi en faisant appel au 2 pour H0 et K0 on peut discriminer les cas pour savoir dans le quel des 3 cas on se trouve
	 LowRankMatrix<T>& hr0 = H0.get_low_rank_block_data();
	 LowRankMatrix<T>& kr0 = K0.get_low_rank_block_data();
	 DenseBlockData<T>& hf0 = H0.get_dense_block_data();
	 DenseBlockData<T>& kf0 = K0.get_dense_block_data();
	//cas facile : les deux sont lowrank on fait direct la multiplication (implémentée?)
	if((hf0 == nullptr) and (kf0 == nullptr)){
	  LowRankMatrix<T> hk0 = h0*k0;}
	//cas 2 : hmat*lowrank
	else if( hr0 == nullptr ){
	  Matrix<T> U = kr0.Get_U(); Matrix<T> V = kr0.Get_V();
	  int rank = U.nb_cols();
	  LowRankMatrix<T> hk0;
	  Matrix<T>* hk0U = hk0.Get_U();
	  Matrix<T>* hk0V = res.Get_V();
	  *hk0V = V;
	  Matrix<T> mattemp(kr0.nb_rows(),U.nb_cols());
	  for ( int k =0; k< rank;++k){
	    vector<T> Uk = U.get_col(k);
	    vector<T> temp = H0*Uk;
	    for ( int l =0; l< U.nb_rows();++l){
	      mattemp(l,k)= Uk[l];
	    }
	  }
	  *hk0U = mattemp;
	}
	//cas 3: on fait le cas deux sur la transposé (et ce qu on peu faire la transposé)
	// Sinon on code un produit vecteur mat xA = At x
	else if( kr0 ==nullptr){
	  
	}
	  
      }
      else{
	pair<Block<T>,Block<T>> HK0; HK0.first = H0;HK0.second = K0;
	Sh0.push_back(HK0);
      }
    }
  }
  std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> S0;
  S0.first = Sr0; S0.second = Sh0;
  return S0;
}  
    MPI_Finalize();
}
