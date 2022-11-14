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

//fonction restrict

// template < typename T>

// std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> restrict(std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>>>> S, Cluster t, Cluster s){
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
//   }
// }
// //Multiplication sans troncature L= HK
// // Block prod(Block,Block) ou hmat prod(hmat,hmat)?
// // quand on appelle la fonction avec t=s=I, S(I,I)= Sh = HxK
// // L'idée c'est de construire L de facon récursive 
// BLock<T> Hmult(Block<T> H, Block<T> K){
//   //On instensie L dans le cas ou H et K sont les matrices entières
//   Cluster H_t = H.get_target_cluster(); Cluster H_s = H.get_source_cluster();
//   Cluster K_t = K.get_target_cluster(); Cluster K_s = K.get_source_cluster();
//   Block<T> H_root = H.get_root() ; Block<T> K_root = K.get_root();
//   Cluster hr_t = H_root.get-target_cluster(); Cluster hr_s = H_root.get_source_cluster();
//   Cluster kr_t = K_root.get_target_cluster() ; Cluster kr_s  = K_root.get_source_cluster();
//   // on test si les cluster sont les mêmes que les racines
//   bool test_H = ((H_t.get_offset() == hr_t.get_offset()) and (H_t.get_size() == hr_s.get size()));
//   bool test_K = ((K_t.get_offset() == kr_t.get_offset()) and (K_s.get_size() == kr_s.get_size()));
//   if ( test_H and test_K){
//     Block<T> L( Admissibility condition, K_t , H_s);
//     //Soit on construit maintenant le block cluster tree soit a chaque étage on fait des build_son
//     L.build('N', false);
//     pair<vector<Block<T>>,pair<vector<Block<T>>,vector<Block<T>>>> S;
//     pair<vector<Block<T>>,vector<Block<T>>>* Sh = &(S.second);
//     pair<vector<Block<T>>,vector<Block<T>>> stemp; stemp.first = H; stemp.second =K; *Sh = stemp; *S.second = Sh;
//   }
//   // Si on est pas tout en haut normalement L et S sont instancié
//   else{
//     // On regarde si tau*sigma = h-t*K_s est feuille de L
//     int repere_L
    


    
// Hmatrix<T> prod(Block<T> L,std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>> S , Cluster t, Cluster s){
//   //Si t,s not a leaf <=> sons(t,s) != emptyset
//   if(!(t.get_sonptr(0)==nullptr) and !(s.get_sonptr(0)==nullptr)){
//     int repere_t = 0;
//     //on regarde les fils de (t,s)
//     while(!(t.get_sonptr(repere_t) == nullptr)){
//       int repere_s = 0;
//       Cluster tson = t.get_sonptr(repere_t);
//       while(!(s.get_sonptr(repere_s) == nullptr)){
// 	Cluster sson = s.get_sonptr(repere_s);
// 	std::pair<std::vector<Block<T>>,std::vector<pair<Block<T>,Block<T>> S0 = restrict(S,tson,sson);
// 	//On 
	
	
	      
    
    
    
  





class MyCondition: public VirtualAdmissibilityCondition {
   bool ComputeAdmissibility(const VirtualCluster &target, const VirtualCluster &source, double eta) const override {
        bool admissible = 2 * std::min(target.get_rad(), source.get_rad()) < eta * std::max((norm2(target.get_ctr() - source.get_ctr()) - target.get_rad() - source.get_rad()), 0.);
        return admissible;
    }
};


// class MyMatrix0 : public VirtualGenerator<double> {
//      const vector<double> &p1;
//     const vector<double> &p2;
//     int space_dim;

//   public:
//     Constructor
//     MyMatrix0(int space_dim0, int nr, int nc, const vector<double> &p10, const vector<double> &p20) : VirtualGenerator(nr, nc), p1(p10), p2(p20), space_dim(space_dim0) {}

//     Virtual function to overload
//     double get_coef(const int &k, const int &j) const {
//       return 0;
//     }

//     Virtual function to overload
//     void copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const override {
//         for (int j = 0; j < M; j++) {
//             for (int k = 0; k < N; k++) {
//                 ptr[j + M * k] = this->get_coef(rows[j], cols[k]);
//             }
//         }
//     }

//     Matrix vector product
//     std::vector<double> operator*(std::vector<double> a) {
//         std::vector<double> result(nr, 0);
//         for (int j = 0; j < nr; j++) {
//             for (int k = 0; k < nc; k++) {
//                 result[j] += this->get_coef(j, k) * a[k];
//             }
//         }
//         return result;
//     }

//     Frobenius norm
//     double norm() {
//         double norm = 0;
//         for (int j = 0; j < nr; j++) {
//             for (int k = 0; k < nc; k++) {
//                 norm += this->get_coef(j, k);
//             }
//         }
//         return norm;
//     }
// };



  ///////////////////////////////

  /////////////////////////////
class MyMatrix : public VirtualGenerator<double> {
    const vector<double> &p1;
    const vector<double> &p2;
    int space_dim;
    vector<double> mat;

  public:
    // Constructor
  MyMatrix(int space_dim0, int nr, int nc, const vector<double> &p10, const vector<double> &p20, const string& matfile) : VirtualGenerator(nr, nc),p1(p10),p2(p20), space_dim(space_dim0) {
    // On load les valeurs de la matrice depuis un .txt//
    ifstream file;
    file.open(matfile);
    vector<double> data;
    //Checking if errors occured
    if(!file){
      std::cerr << "Error : Can't open the file" << std::endl;}
    else{
      //Creating a word string
      std::string word;
      //Creating a vector of string for file content
      std::vector<std::string> file_content;
         
      int count = 0;
         
      while(file >> word)
      {
      //Adding every words to the file content
	file_content.push_back(word);
      //Increasing the word counter;
	count++;}
      int k=0 ; int l=0;
      //Printing file content
        for(int i = 0; i < count; i++){
	  string str = file_content[i];	  double d1;
	  stringstream stream1;
	  stream1 << str;
	  stream1 >> d1;
	  data.push_back(d1);
        }
    }  
    //Closing the file
    file.close();
    mat = data;
  }


  //   Virtual function to overload
    double get_coef(const int &k, const int &j) const {
       return mat[nr*k+j];
    }

  // Virtual function to overload
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const override {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                ptr[j + M * k] = this->get_coef(rows[j], cols[k]);
            }
        }
    }

  //Matrix vector product
    std::vector<double> operator*(std::vector<double> a) {
        std::vector<double> result(nr, 0);
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                result[j] += this->get_coef(j, k) * a[k];
            }
        }
        return result;
    }

  //Frobenius norm
    double norm() {
        double norm = 0;
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
	      norm += pow(this->get_coef(j, k),2);
            }
        }
        return sqrt(norm);
    }


    // pour load la matrice de xavier
   //  ifstream file;
  //   file.open(matfile);
  //   vector<complex<double>> data;
  //   //Checking if errors occured
  //   if(!file){
  //     std::cerr << "Error : Can't open the file" << std::endl;}
  //   else{
  //     //Creating a word string
  //     std::string word;
  //     //Creating a vector of string for file content
  //     std::vector<std::string> file_content;
         
  //     int count = 0;
         
  //     while(file >> word)
  //     {
  //     //Adding every words to the file content
  // 	file_content.push_back(word);
  //     //Increasing the word counter;
  // 	count++;}
  //     int k=0 ; int l=0;
  //     //Printing file content
  //       for(int i = 0; i < count; i++){
  // 	  string str = file_content[i];
  // 	  double d;
  // 	  string str1; string str2;
  // 	  int kk = 1;
  // 	  while(str[kk] != ','){
  // 	    str1.push_back(str[kk]);
  // 	    kk+=1;}
  // 	  kk+=1;
  // 	  double d1;
  // 	  stringstream stream1;
  // 	  stream1 << str1;
  // 	  stream1>> d1;
  // 	  while(str[kk] != ')'){
  // 	    str2.push_back(str[kk]);
  // 	    kk+=1;}
  // 	  double d2;
  // 	  stringstream stream2;
  // 	  stream2<< str2;
  // 	  stream2 >> d2;
  // 	  complex<double> z (d1,d2);
  // 	  data.push_back(z);
  //       }
  //   }  
  //   //Closing the file
  //   file.close();
  //   mat = data;
  // }

//     // Virtual function to overload
//   complex<double> get_coef(const int &k, const int &j) const {return mat[nr*k+j];}

//     // Virtual function to overload
// void copy_submatrix(int M, int N, const int *const rows, const int *const cols, complex<double> *ptr) const override {
//         for (int j = 0; j < M; j++) {
//             for (int k = 0; k < N; k++) {
//                 ptr[j + M * k] = this->get_coef(rows[j], cols[k]);
//             }
//         }
//     }

//     // Matrix vector product
// std::vector<double> operator*(std::vector<complex<double>> a) {
//   std::vector<double> result(nr, 0);
//         for (int j = 0; j < nr; j++) {
//             for (int k = 0; k < nc; k++) {
// 	      complex<double> z = this->get_coef(j,k);
// 	      complex<double> aa = a[k];
// 	      double zx = z.real();double zy = z.imag(); double ax = aa.real(); double ay = aa.imag();
//               result[j] += z*aa;
//             }
//         }
//         return result;
//     }

//   std::vector<complex<double>> operator*(std::vector<double> a) {
//   std::vector<complex<double>> result(nr, 0);
//         for (int j = 0; j < nr; j++) {
//             for (int k = 0; k < nc; k++) {
// 	      complex<double> z = this->get_coef(j,k);
// 	      double zx = z.real();double zy = z.imag();
//               result[j] += z*a[k];
//             }
//         }
//         return result;
//     }

//     // Frobenius norm
//     double norme() {
//         double n = 0;
//         for (int j = 0; j < nr; j++) {
//             for (int k = 0; k < nc; k++) {
// 	      complex<double> z = this->get_coef(j,k);
// 	      n += norm(z);
//             }
//         }
//         return n;
//     }
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

    //std::string outputpath = argv[1];

    // Htool parameters
    double epsilon = 0.0000001;
    double eta     = 1;

    // n² points on a regular grid in a square
    int n    = 6021;
    int size = 15030;

    // // p1: points in a square in the plane z=z1
    // double z = 1;
    // vector<double> p(3 * size);
    // for (int j = 0; j < n; j++) {
    //     for (int k = 0; k < n; k++) {
    //         p[3 * (j + k * n) + 0] = j;
    //         p[3 * (j + k * n) + 1] = k;
    //         p[3 * (j + k * n) + 2] = z;
    //     }
    // }
        // // p1: points su z=z1
    double z = 1;
    vector<double> p(3 * size);
    for (int j = 0; j < size; j++) {
            p[3 * j + 0] = j;
            p[3 * j + 1] = z;
            p[3 *j + 2] = z;
        }
  // vector<double> nodes;
  // vector<double> Nodes;
  // std::ifstream file;
  // file.open("/work/sauniear/Documents/matrice_Xavier/maillage_arthur.txt");
  //   //Checking if errors occured
  // if(!file){
  //   std::cerr << "Error : Can't open the file" << std::endl;}
  //   else{
  //     //Creating a word string
  //     std::string word;
  //     //Creating a vector of string for file content
  //     std::vector<std::string> file_content;  
  //     int count = 0;   
  //     while(file >> word)
  //     {
  //     //Adding every words to the file content
  // 	file_content.push_back(word);
  //     //Increasing the word counter;
  // 	count++;}
  //     //Printing file content
  //       for(int i = 0; i < count; i++){
  // 	  string str = file_content[i];
  // 	  double p;
  // 	  stringstream stream;
  // 	  stream << str;
  // 	  stream >> p;
  // 	  nodes.push_back(p);
  //       }
  // 	 for (int k = 0; k<6021*6021; ++k){
  //   double d1 = nodes[4*k+1+7]; double d2 = nodes[4*k+2+7]; double d3 = nodes[4*k+3+7];
  //   Nodes.push_back(d1); Nodes.push_back(d2); Nodes.push_back(d3);}
  //       }
  //   //Closing the file}
  // file.close();
  // vector<double> p = Nodes;


  // vector<double> nodes;
  // vector<double> Nodes;
  // std::ifstream file;
  // file.open("/work/sauniear/Documents/matrice_Xavier/maillage_arthur.txt");
  //   //Checking if errors occured
  // if(!file){
  //   std::cerr << "Error : Can't open the file" << std::endl;}
  //   else{
  //     //Creating a word string
  //     std::string word;
  //     //Creating a vector of string for file content
  //     std::vector<std::string> file_content;  
  //     int count = 0;   
  //     while(file >> word)
  //     {
  //     //Adding every words to the file content
  // 	file_content.push_back(word);
  //     //Increasing the word counter;
  // 	count++;} 
  //     //Printing file content
  //       for(int i = 0; i < count; i++){
  // 	  string str = file_content[i];
  // 	  double p;
  // 	  stringstream stream;
  // 	  stream << str;
  // 	  stream >> p;
  // 	  nodes.push_back(p);
  //       }
  // 	int nn = (nodes.size()-7)/4;
  // 	 for (int k = 0; k<nn; ++k){
  //   double d1 = nodes[4*k+1+7]; double d2 = nodes[4*k+2+7]; double d3 = nodes[4*k+3+7];
  //   Nodes.push_back(d1); Nodes.push_back(d2); Nodes.push_back(d3);}
  //       }
  // //   //Closing the file}
  // file.close();
  // Nodes.resize(3*6021);
  // vector<double> p = Nodes;
    // Hmatrix
MyMatrix A(3,size,size,p,p,"/work/sauniear/Documents/matrice_Abdoulaye/Hmatrix/test_F3.txt");
std::vector<double> x(size, 1);
std::vector<double> result(size, 0);
std::shared_ptr<Cluster<PCA<SplittingTypes::RegularSplitting>>> t = make_shared<Cluster<PCA<SplittingTypes::RegularSplitting>>>(3);
t->build(size,p.data(), 2);
HMatrix<double> HA(t,t,epsilon, eta, 'N', 'N');
HA.build(A, p.data());
 result = HA * x;

    // Output
    HA.print_infos();
    HA.save_plot("mat_art_plot");
    // HA.get_target_cluster()->save_geometry(p.data(), outputpath + "/smallest_example_cluster", {1, 2, 3});
    // std::cout << outputpath + "/smallest_example_plot" << std::endl;
    std::cout << Frobenius_absolute_error(HA, A) / A.norm() << std::endl;
    std::cout << norm2(A * x - result) / norm2(A * x) << std::endl;
    std::cout<< A.get_coef(1,1)<<','<<A.get_coef(3,8)<<endl;





























    // std::vector<std::pair<int,int>> ix = HA.get_MasterOffset_t();
    // std::vector<std::pair<int,int>> iy = HA.get_MasterOffset_s();
    // cout<< "%%%%%%%%%%%%%%%%%taille%%%%%%%%%%%"<<endl;
    // cout<< ix.size() << endl;
    // cout<< iy.size() << endl;
    // for (int k =0 ; k<10; ++k){
    //   pair<int,int> ixk = ix[k]; pair<int,int> iyk = iy[k];
    //   cout<< "t:"<< ixk.first<< ","<< ixk.second<<"%%%%%%%"<< "s:"<< iyk.first<<","<<iyk.second<< endl;}

    //   std::vector<SubMatrix<std::complex<double>> *> nf = HA.get_MyNearFieldMats();
    //   std::vector<LowRankMatrix<std::complex<double>> *> ff = HA.get_MyFarFieldMats();
    //   std::vector<SubMatrix<std::complex<double>> *> dnf = HA.get_MyDiagNearFieldMats();
    //   std::vector<LowRankMatrix<std::complex<double>> *> dff = HA.get_MyDiagFarFieldMats();
    //   std::vector<SubMatrix<std::complex<double>> *> sdnf = HA.get_MyStrictlyDiagNearFieldMats();
    //   std::vector<LowRankMatrix<std::complex<double>> *> sdff = HA.get_MyStrictlyDiagFarFieldMats();

    //   cout<< "$$$$$$$$$$$$$"<<"test sur les get"<<"$$$$$$$$$$$"<< endl;
    //   cout<< nf.size()<< endl;
    //   cout<< ff.size()<< endl;
    //   cout<< dnf.size() << endl;
    //   cout<< dff.size()<< endl;
    //   cout<< sdnf.size()<< endl;
    //   cout<< sdff.size()<< endl;

      // cout << "test sur les blocs"<< endl;
      // for (int k =0; k <10;++k){
      // 	LowRankMatrix<std::complex<double>> lr =*ff[k];
      // 	cout<< "étape"<< k<< endl;
      // 	cout<< "rank"<< lr.rank_of() << endl;
      // 	cout <<"get_ir.size"<< (lr.get_ir()).size()<< endl;
      // 	cout<<"get_ir"<< (lr.get_ir())[0]<<","<<(lr.get_ir())[1]<<","<<(lr.get_ir())[2]<< endl;
      // 	cout <<"get_ic.size"<< (lr.get_ic()).size()<< endl;
      // 	cout<<"get_ic"<< (lr.get_ic())[0]<<","<<(lr.get_ic())[1]<<","<<(lr.get_ic())[2]<< endl;}

      // cout<< "test sur la hmatrice"<< endl;
      // vector<pair<int,int>> mt = HA.get_MasterOffset_t();
      // vector<pair<int,int>> ms = HA.get_MasterOffset_s();
      // int nnn = mt.size();
      // cout<< "nb d'offset " << nnn<< endl;
      // for (int k =0; k<10; ++k){
      // 	pair<int,int> pt, ps; pt = mt[k];ps=ms[k];
      // 	cout<< "offset "<<k<<" t:"<< pt.first<<","<<pt.second<<" ,s:"<<ps.first<<","<<ps.second << endl;
      // }
      // cout<< "test des offset pour les far fields"<< endl;
      
      // for (int kkk =0; kkk< 10;++kkk){
      // 	LowRankMatrix<complex<double>> tempo = *ff[kkk];
      // 	cout<<"taille:"<< (tempo.get_ir()).size()<<','<< (tempo.get_ic()).size()<< endl;
      // 	cout<<"i:"<< tempo.get_offset_i()<<",j:"<< tempo.get_offset_j()<< endl;}

      // cout<< "test des offset pour les near fields"<< endl;
      // for(int kkkk =0; kkkk<10;++kkkk){
      // 	SubMatrix<complex<double>> tempo = *nf[kkkk];
      // 	cout<<"taille:"<< (tempo.get_ir()).size()<<','<< (tempo.get_ic()).size()<< endl;
      // 	cout<<"i:"<< tempo.get_offset_i()<<" ,j:"<< tempo.get_offset_j()<< endl;
      // }


      // cout<< "test des offset pour les diag far fields"<< endl;
      
      // for (int k =0; k< 10;++k){
      // 	LowRankMatrix<complex<double>> tempo = *dff[k];
      // 	cout<<"i:"<< tempo.get_offset_i()<<",j:"<< tempo.get_offset_j()<< endl;}

      // cout<< "test des offset pour les diag near fields"<< endl;
      // for(int k =0; k<10;++k){
      // 	SubMatrix<complex<double>> tempo = *dnf[k];
      // 	cout<<"i:"<< tempo.get_offset_i()<<" ,j:"<< tempo.get_offset_j()<< endl;
      // }


      // cout<< "test des offset pour les strictly diag near fields"<< endl;
      // for(int k =0; k<10;++k){
      // 	SubMatrix<complex<double>> tempo = *sdnf[k];
      // 	cout<<"i:"<< tempo.get_offset_i()<<" ,j:"<< tempo.get_offset_j()<< endl;
      // }


	

      // 	LowRankMatrix<std::complex<double>>* lr =ff[0];
      // 	cout<< "test low rank"<< endl;
      // 	cout<< lr->get_U(0,0)<< endl;complex<double> z (3.14,3.14);
      // 	complex<double>* zz = &z;
      // 	lr->assign_U(0,0,zz);
      // 	LowRankMatrix<std::complex<double>> llr =*ff[0];
      // 	cout << llr.get_U(0,0)<< endl;
      // cout<<"teeeeeeeeeeeeeest"<<endl;
      // LowRankMatrix<std::complex<double>>* lr0 = ff[0];
      // LowRankMatrix<std::complex<double>> lr1 = *ff[1];
      // cout<< lr0->get_U(0,0)<<"'"<< lr0->rank_of()<<endl;
      // *lr0 = lr1;
      // cout<<  lr0->get_U(0,0)<<"'"<< lr0->rank_of()<<endl;
      // std::vector<LowRankMatrix<std::complex<double>> *> fff = HA.get_MyFarFieldMats();
      //  LowRankMatrix<std::complex<double>>* lr00 = fff[0];
      //  cout<<  lr0->get_U(0,0)<<"'"<< lr0->rank_of()<<endl;


      //  cout<< "teeeeeeeeeeest pour les full"<< endl;
       // SubMatrix<complex<double>>* sr0 = nf[0];
       // // SubMatrix<complex<double>> temp(*sr0);
       // cout<< (*sr0)(0,0)<< endl;
       // complex<double> z (3.14,3.14);
       // (*sr0)(0,0) = z;
       // std::vector<SubMatrix<std::complex<double>> *> sss = HA.get_MyNearFieldMats();
       // SubMatrix<complex<double>>* sr00 = sss[0];
       // cout<< (*sr0)(0,0)<<endl;
       // cout<< (*sr00)(0,0)<< endl;
       
      
      
	


    


    
    	// LowRankMatrix<std::complex<double>> teest =*ff[0];
	// vector<int> irr = teest.get_ir();
	// cout<< teest.get_offset_i()<< endl;
	// for (auto it = irr.begin(); it != irr.end(); ++it){
	//   cout<< *it<< endl;}
	// cout<<"-------------------------"<< endl;
	// cout<< teest.get_offset_j()<< endl;
	// vector<int> icc = teest.get_ic();
	// for (auto it = icc.begin(); it != icc.end(); ++it){
	//   cout<< *it<< endl;}
	// cout<<"***********************"<< endl;
	// cout<<"xr"<< endl;
	// vector<int> xxr; vector< int> xxc;
	// xxr = teest.get_xr(); xxc = teest.get_xc();
	// for ( auto it = xxr.begin(); it != xxr.end();++it){
	//   cout<< *it<< endl;}
	// cout<< "xc"<< endl;
	// for (auto it = xxc.begin(); it != xxc.end; ++it){
	//   cout<< *it << endl;}
	// cout<<*teest<< endl;
	// int r = teest.rank_of();
        // int nr = irr.size(); int nc = icc.size();
	// Matrix<complex<double>> UU(nr, r);
	// Matrix<complex<double>> VV(r,nc);
	// for (int k =0; k< nr;++k){
	//   for ( int l = 0; l< r ; ++l){
	//     complex<double> z = teest.get_U(k,l);
	//     UU(k,l) = z;
	//   }
	// }
	// for(int k =0; k< r;++k){
	//   for (int l =0;l<nc; ++l){
	//     complex<double> z = teest.get_V(k,l);
	//     VV(k,l) = z;}
	// }
        // Matrix<complex<double>> prod = UU*VV;
	// for (int k=0; k < nr; ++k){
	//   for(int l =0; l<nc;++l){
	//     cout<< prod(k,l)<<'\t';}
	//   cout<< '\n';}

    // 	cout<< "test block"<< endl;
    //     Cluster<PCA<SplittingTypes::RegularSplitting>> temp = Cluster<PCA<SplittingTypes::RegularSplitting>>(3);
    // temp.build(size,p.data(),2);
    // MyCondition *m;
    // RjasanowSteinbach* mm;
    // Block B(mm,temp,temp);
	// B.build_son(t,t);
	// Block b = B.get_son(0);
	// for (auto it = t.begin(); it != t.end();++it) {
	//   cout <<*it<< endl;}
        // std::shared_ptr<Cluster<PCA<SplittingTypes::RegularSplitting>>> tt = b.get_target_cluster();
	// for (auto it = tt.begin(); tt!= tt.end(); ++it){
	//   cout<< *it<< endl;}
    // B.build_son(temp,temp);
    // for (int k =0; k<1 ; ++k){
    //   Block &bb = B.get_son(k);
    //   auto tt = bb.get_target_cluster();
    //   // Cluster<PCA<SplittingTypes::RegularSplitting>> ty = bb.get_cluster_source();
    //   // Cluster tx = bb.get_target_cluster();
    //   // cout << bb.get_target_cluster()<< endl;
    //   // cout << bb.get_source_cluster()<< endl;
    //   cout<< "************************"<< endl;}
    // B.build('N', false, MPI_COMM_WORLD);
    // const std::vector<Block*> vb = B.get_tasks();
    // cout<< vb.size()<< endl;
    // B.build_son(B.get_target_cluster() ,B.get_source_cluster() );
    // cout<< B.get_eta() << endl;
    // cout<< B.IsAdmissible() << endl;
    // Cluster<PCA<SplittingTypes:: RegularSplitting>> tt ( B.get_target_cluster());
    // Cluster<PCA<SplittingTypes::RegularSplitting>> ss (  B.get_source_cluster());
    // auto t = B.get_source_cluster(); auto s =B.get_target_cluster();
    // B.build_son(tt,ss);
    // B.build('N');

    // vector<Block*> vb = B.get_tasks();
    // Finalize the MPI environment.
    // cout<< "-----------------"<< endl;
    // cout<< "test cluster"<< endl;
    // cout<< "----------------"<< endl;
    // Ca marche et oin peut pointer vers les sons il reste a lier les cluster a leur blocs
    
    // int taille = t->get_nb_sons();
    // cout<<"étape"<<' '<<0<< endl;
    // cout<< "nb sons "<<' '<< taille<< endl;
    // cout<< "depth "<<t->get_depth()<< endl;
    // cout<< "rank "<<t->get_rank() << endl;
    // cout<< "offset"<<t->get_offset() << endl;
    // cout<< "centre "<<t->get_ctr()<<endl;
    // for (int k=0; k< taille; ++k){
    //   cout<< "étape"<< ' '<<k+1<< endl;
    //   auto a =t->get_son_ptr(k);
    //   cout<<"nb sons "<< a->get_nb_sons()<<endl;
    //   cout<<"depth "<< a->get_depth()<< endl;
    //   cout<<"rank "<<  a->get_rank() << endl;
    //   cout<< "offset "<<a->get_offset() << endl;
    //   cout<< "centre "<<a->get_ctr()<<endl;
    //   for (int l=0;l< a->get_nb_sons();++l){
    // 	cout<< "étape"<< ' '<<k+1<<' '<<l<< endl;
    // 	 auto aa =a->get_son_ptr(l);
    // 	 cout<<"nb sons "<< aa->get_nb_sons()<<endl;
    // 	 cout<<"depth "<< aa->get_depth()<< endl;
    // 	 cout<<"rank "<<  aa->get_rank() << endl;
    // 	 cout<< "offset "<<aa->get_offset() << endl;
    // 	 cout<< "centre "<<aa->get_ctr()<<endl;
    //   }}

    //Test du truc que j'ai implémenter pour avoir le blocktree

    // std::unique_ptr<Block<complex<double>>> test_block = HA.get_BlockTree();
    // std::vector<Block<complex<double>>*> test_computed = HA.get_ComputedBlock();
    // cout<< test_computed.size()<< endl;
    // for (int k =0;k<20; ++k){
    //   Block<complex<double>>* B1 = test_computed[k];
    //   // VirtualCluster v1 = B1->get_target_cluster();
    //   // VirtualCluster v2 = B1->get_source_cluster();
    //   // int a1 = B1-> get_mintargetdepth();
    //   // int a2 = B1-> get_minsourcedepth();
    //   // cout<< a1<< "et"<<a2<< endl;
    //   cout<<"taille"<<' '<<B1->get_size()<< endl;
    //   cout<<"offset des blocs"<< endl;
    //   Block<complex<double>> B2 = *B1;
    //   // int o1 = B1->get_target_cluster().get_offset();
    //   // int o2 = B1->get_source_cluster().get_offset();
    //   // cout<< o1<< endl;
    //   // cout<< o2<< endl;
    //         int o1 = B1->get_target_cluster().get_offset();
    //   int o2 = B1->get_source_cluster().get_offset();
    //   cout<< o1<< endl;
    //   cout<< o2<< endl;
    //   cout<<"-------------"<<endl;
    //   const DenseBlockData<complex<double>>* db = B1->get_dense_block_data();
    //   const LowRankMatrix<complex<double>>* lr = B1->get_low_rank_block_data();
    //   if (db == nullptr){
    // 	// cout<<"xr: ";
    // 	// vector<int> xr = lr->get_xr();
    // 	// for (auto it = xr.begin(); it!=xr.end();++it){
    // 	//   cout<< *it<< "\t";
    // 	// }
    // 	// cout<<endl;
    // 	// cout<<"xc: ";
    // 	// vector<int> xc = lr->get_xc();
    // 	// for (auto it = xc.begin(); xc!=xc.end();++it){
    // 	//   cout<<xc[k]<<"\t";}
    // 	// cout << endl;}
    // 	int n1= lr->nb_rows(); int n2 = lr->nb_cols();
    // 	cout<<"nr,nc= "<< n1<<','<<n2<< endl;
    // 	// const Matrix<complex<double>> MM(lr->nb_rows(),lr->nb_cols());
    //     // complex<double> z (0,0);
    //     // complex<double>* zz =&z;
    // 	// lr->get_whole_matrix(zz);
    // 	// cout<<*zz<< endl;
    // 	// cout<<"------------------------"<<endl;
    //   }
    //   // else{
    //   // 	// cout<<"get_row"<<endl;
    //   // 	// vector<int> xx = db->get_row();
    //   // 	// vector<int> yy =db->get_col();
    //   // 	// for ( auto it = xx.begin();it!=xx.end();++it){
    //   // 	//   cout<< *it<<"\t";}
    //   // 	// cout<< endl;
    //   // 	// cout<<"get_col"<< endl;
    //   // 	// for (auto it = yy.begin(); it!= yy.end();++it){
    //   // 	//   cout<<*it<<"\t";
    //   // 	// }
    //   // 	// cout<< endl;
    //   // // 	cout<<"print matrice";
    //   // // 	for(int k =0; k< db->nb_rows();++k){
    //   // // 	  for (int l =0; l< db->nb_cols();++l){
    //   // // 	    cout<< (*db)(k,l)<<"\t";}
    //   // // 	  cout<<endl;}
    //   // // 	cout<<"-------------------------"<< endl;
    //   // }
    //   }
    // string Myfile("/work/sauniear/Documents/test_C++/offset_H.txt");
    // ofstream flux(Myfile.c_str());
    // for (int k =0;k< test_computed.size();++k){
    //   Block<complex<double>>* B1 = test_computed[k];
    //   cout<<"taille"<<' '<<B1->get_size()<< endl;
    //   cout<<"offset des blocs"<<"\t";
    //   int o1 = B1->get_target_cluster().get_offset();
    //   int o2 = B1->get_source_cluster().get_offset();
    //   cout<< o1<< endl;
    //   cout<< o2<< endl;
    //   int nn1 = B1->get_target_cluster().get_size();
    //   int nn2 = B1->get_source_cluster().get_size();
    //   flux<<nn1;flux<<' ';flux<<o1; flux<<' '; flux<< nn2; flux<<' ';flux<< o2; flux<< endl;}

    // Test pour descendre dans la matrice

    // Block <complex<double>>* B3 = test_computed[0];
    // Block <complex<double>>* rB = B3->get_root();
    // int a1 = B3->get_target_cluster().get_size();
    // int a2 = B3->get_source_cluster().get_size();
    // int b1 = rB->get_target_cluster().get_size();
    // int b2 = rB->get_source_cluster().get_size();
    // cout<< "test pour la root"<<endl;
    // cout<<"fille"<<' '<<a1<<','<<a2<<endl;
    // cout<<"mère"<<' '<<b1<<','<<b2<<endl;


    // Block<complex<double>> &fille1 = rB->get_son(0);
    // int c1 = fille1.get_target_cluster().get_offset();
    // int c2 = fille1.get_source_cluster().get_offset();
    // cout<<"étage 1"<<' '<<c1<<','<<c2<<endl;
    // Block<complex<double>> &fille2 = fille1.get_son(0);
    // int c11 = fille2.get_target_cluster().get_offset();
    // int c21 = fille2.get_source_cluster().get_offset();
    // cout<<"étage 2"<<' '<<c11<<','<<c21<<endl;
    // Block<complex<double>> &fille3 = fille1.get_son(1);
    // int c12 = fille3.get_target_cluster().get_offset();
    // int c22 = fille3.get_source_cluster().get_offset();
    // cout<<"étage 2"<<' '<<c12<<','<<c22<<endl;
    // Block<complex<double>> &fille4 = fille1.get_son(2);
    // int c13 = fille4.get_target_cluster().get_offset();
    // int c23 = fille4.get_source_cluster().get_offset();
    // cout<<"étage 2"<<' '<<c13<<','<<c23<<endl;
    // Block<complex<double>> &fille5 = fille1.get_son(3);
    // int c14 = fille5.get_target_cluster().get_offset();
    // int c24 = fille5.get_source_cluster().get_offset();
    // cout<<"étage 2"<<' '<<c14<<','<<c24<<endl;
    // int rep = 0;
    // while( !(&fille1.get_son(rep) ==  nullptr)){
    //   Block<complex<double>> & fille = fille1.get_son(rep);
    //   int rep1 = fille.get_target_cluster().get_size();
    //   int rep2 = fille.get_source_cluster().get_size();
    //   int ofs1 = fille.get_target_cluster().get_offset();
    //   int ofs2 = fille.get_source_cluster().get_offset();
    //   cout<<"étage 2"<<endl;
    //   cout<<"taille"<<rep1<<','<<rep2<<endl;
    //   cout<<"offset"<<ofs1<<','<<ofs2<< endl;
    //   int rep0 = 0;
    //   while(!(&fille.get_son(rep0) == nullptr)){
    // 	 Block<complex<double>> & fille0 = fille.get_son(rep);
    // 	 int rep10 = fille0.get_target_cluster().get_size();
    // 	 int rep20 = fille0.get_source_cluster().get_size();
    // 	 int ofs10 = fille0.get_target_cluster().get_offset();
    // 	 int ofs20 = fille0.get_source_cluster().get_offset();
    // 	 cout<<"étage 3"<<endl;
    // 	 cout<<"taille"<<rep10<<','<<rep20<<endl;
    // 	 cout<<"offset"<<ofs10<<','<<ofs20<< endl;
    // 	 rep0 +=1;}
    //   rep+= 1;}
      
	  






    // test pour la version structural changes




    // test pour les offset



    // cout<< "test pour les offset"<< endl;
    // cout<<"********"<< endl;
    // std::vector<pair<int,int>> oft = HA.get_MasterOffset_t();
    // std::vector<pair<int,int>> ofs = HA.get_MasterOffset_s();
    // for(int k =0; k<20 ; ++k){
    //   pair<int,int> oftt = oft[k]; pair<int,int> ofss = ofs[k];
    //   cout<< "t="<< oftt.first<<","<<oftt.second<< " ,s="<<ofss.first<< ","<< ofss.second<< endl;
    // }


    // // Matrix<complex<double>> M1 = HA.get_local_interaction();
    // Matrix<complex<double>> M2 = HA.get_local_diagonal_block();
    // // cout<< normFrob(M1-M2)<<endl;
    // // cout<<"size"<< endl;
    // // cout<< M1.nb_cols()<<','<< M1.nb_rows()<<" ,"<<M2.nb_cols()<<M2.nb_rows()<< endl;
    // cout<< M2.nb_rows()<<','<< M2.nb_cols()<< endl;


    //visualisation des offset

    // string Myfile("/work/sauniear/Documents/test_C++/offset_H.txt");
    // ofstream flux(Myfile.c_str());
    // cout<< ofs.size()<< endl;
    // cout<< sizeof(ofs) << endl;
    // for (int k =0 ;k<10000000;++k){
    //   pair<int,int> oftt = oft[k]; pair<int,int> ofss = ofs[k];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;}


    // int k1 =0;
    // for (auto it = oft.begin(); it!=oft.end();++it){
    //   pair<int,int> oftt = oft[k1]; pair<int,int> ofss = ofs[k1];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;
    //   k1=k1+1;}
    // cout<<k1<< endl;
    // while(oft[k1] != ){
    //         pair<int,int> oftt = oft[k1]; pair<int,int> ofss = ofs[k1];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;
    //   k1=k1+1;}
    // cout <<oft::size<< endl;
    // for (int k1=0; k1<oft::size;++k1){
    //   pair<int,int> oftt = oft[k1]; pair<int,int> ofss = ofs[k1];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;}
    // int k1=0;
    // for (auto it=oft.begin(); it!=oft.end();++it){
    //   pair<int,int> oftt = *it; pair<int,int> ofss = ofs[k1];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;
    //   k1=k1+1;}
    // int k1 = 0;
    // while( oft[k1] != (oft.data).end()){
    //         pair<int,int> oftt = oft[k1]; pair<int,int> ofss = ofs[k1];
    //   flux<< oftt.first<<' '<<oftt.second<<' '<<ofss.first<<' '<<ofss.second<< endl;
    //   k1=k1+1;}
    // cout<< k1<< endl;



    // cout<< k1<<endl;

    cout<< "test pour le produit de lrmat"<< endl;
    // std::vector<LowRankMatrix<std::complex<double>> *> ff = HA.get_MyFarFieldMats();
    // lowRankMatrix<std::complex<double>> M1= *ff[1];
    // lowRankMatrix<std::complex<double>> M2= *ff[2];
    // lowRankMatrix<std::complex<double>> M3 =M1;
    // Block<complex<double>> A1,B1,A2,B2;
    // Block<complex<double>>* A3,B3
    //   A1= M1.Get_U(); A2=M2.get_U(); A3 = &M3.Get_U();
    // B1= M1.Get_V(); B2= M2.Get_V(); B3 = &M3.get_V();
    // *A3 = A1*A2;
    // *B3 = B1*B2;
    // cout<< A1<< endl;
    // cout<< A3 << endl;
    // On va essayer de modifier de la matrice L = HK initialiser a 0


    
    // MyMatrix0 MM(3,size, size,p,p);

    // HMatrix<complex<double>> Hres(t, t, epsilon, eta, 'N', 'N');
    // Hres.build(A, p.data());
    // vector<Block<complex<double>>*> CB = HA.get_ComputedBlock();
    // vector<Block<complex<double>>*> lr;
    // vector<Block<complex<double>>*> full;
    // for( auto it = CB.begin(); it!= CB.end(); ++it){
    //   Block<complex<double>>* bb = *it;
    //   int rep =  bb->get_rank_of();
    //   if( rep == -1){
    // 	full.push_back(bb);}
    //   else{
    // 	lr.push_back(bb);}
    // }


    // Block<complex<double>>* m1= lr[45]; Block<complex<double>>* m2 = lr[1];
    // cout<< m1->get_target_cluster().get_offset()<<endl;
    // Block<complex<double>> mm1 = *m1 ;Block<complex<double>> mm2 = *mm2;
    // LowRankMatrix<complex<double>> a1 = *(m1->get_low_rank_block_data());
    // LowRankMatrix<complex<double>> a2 =*(m2->get_low_rank_block_data());
    // Cluster t_res = *(m1->get_target_cluster()); Cluster  s_res  = *(m2->get_source_clster());
    // int of_t = t_res.get_offset(); int  sz_t = t_res.get_sise();
    // int of_s = ss_res.get_size(); it sz_s = s_temp.get_size();
    // vector<Block<complex<double>>*> b_res = MM.ComputedBlock();
    // Block<complex<double>> b_temp
    // for (auto it= b_res.begin(); it != b_res.end(); ++it){
    //   Cluster t_temp =  it->get_terget_cluster(); Cluster s_temp = it-> get_source_cluster();
    //   int of_temp = t_temp.get_offset(); int  sz_temp = t_temp.res.get_sise();
    //   int of_semp = s_temp.get_size(); it sz_semp = s_temp.get_size();
    //   if ( (of_temp == of_t) and ( ssz_s ==sz_temp)){
    // 	btemp = *it;
    // 	break;
    //   }
    // }
    // LowRankMatrix<complex<double>> lrtarg = btemp.get_low_rank_block_data();
    // Matrix<complex<double>>& tragU = lrtarg.Get_U();Matrix<complex<double>>& tragV = lrtarg.Get_V();
    // *targU = a1 ; *targV = b1*(a2*b3);


    //
      
      
    
										      
    
      
      
    


    
    MPI_Finalize();
    return 0;
}
