#include <htool/htool.hpp>
#include <memory>
#include <iomanip>
using namespace std;
using namespace htool;


// Condition d'adm 
class Mycondition final : public VirtualAdmissibilityCondition {
  public:
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

//Fonction pour découper les blocs non admissibles

void Split (Block<double>& L){
	int rep = 0;
	for (int k = 0 ; k < L.get_target_cluster().get_nb_sons(); ++k){
		for ( int l=0 ; l< L.get_source_cluster().get_nb_sons(); ++l){
			L.build_son(L.get_target_cluster().get_son(k), L.get_source_cluster().get_son(l));
			if (!L.get_son(rep).IsAdmissible()){
				Split(L.get_son(rep));}
			rep +=1;} } 
}
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



    // fonction restrict : les offsets sont bons, il faut implémenter eveluate pour être sure qu'elle marche
    // Au début ca marche bien mais a la fin il me met des blocs pas admissible ( style ofset_t = ofset_s = 0) du coup je comprend pas trop
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
	    //cout<< "alors c'est comment:"<< U.nb_rows()<<','<< U.nb_cols()<<','<<V.nb_rows()<<','<<V.nb_cols()<< endl;
	    //cout<< oft<<','<<ofs<<','<<ofu<<','<<ofv<<endl;
	    //cout<<szs<<','<<szt<<endl;
            if( (oft >=ofu) and ( szt+(oft-ofu) < U.nb_rows()) and (ofs >= ofv) and ( szs+(ofs-ofv) < V.nb_cols())){
                    // si on  a des blocs concernés on extrait
                    // ------> implémenter un extract rapide
                    Matrix<double> Uk(szt, U.nb_cols()); Matrix<double> Vk(V.nb_rows(), szs);
                    for (int i = 0 ; i< szt; ++i){
                        for(int j =0; j < U.nb_cols(); ++j){
                            Uk(i,j) = U(oft-ofu+i,j);
                        }
                    }
                    for (int i = 0 ; i< V.nb_rows(); ++i){
                        for(int j =0; j < szs; ++j){
                            Vk(i,j) = V(i, ofs-ofv+j);
                        }
                    }
                    sr.push_back(Uk);sr.push_back(Vk); of.push_back(oft); of.push_back(ofs);}
                }

	     cout<< "aa"<< endl;

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
                    // cout<< H->nb_sons()<<","<<K->nb_sons()<<endl;
                    // for (int i = 0; i< H->nb_sons();++i){
                    //     Block<double>& Hk = H->get_son(i);
                    //     for(int j=0; j< K->nb_sons(); ++j){
                    //         Block<double>& Kk = K->get_son(j);
		    // 	    if ( (oft == Hk.get_target_cluster().get_offset()) and (ofs == Kk.get_source_cluster().get_offset()) and (szt == Hk.get_target_cluster().get_size()) and (szs== Kk.get_source_cluster().get_size()) ){
		    // 		/*	if ((Hk.get_source_cluster().get_offset() == Kk.get_targer 
		    // 			if ( ((oft>=Hk.get_target_cluster().get_offset() ) and (szt<= Hk.get_target_cluster().get_size())) and ((ofs>=Kk.get_source_cluster().get_offset() ) and (szs<= Kk.get_source_cluster().get_size()))){ */
		    // 	      cout<<'!'<<i<<','<<j<<endl;
                    //             if( (Hk.get_source_cluster().get_size()==Kk.get_target_cluster().get_size())){
                    //                 vtemp.push_back(&Hk); vtemp.push_back(&Kk); }	} }

                    // }} }
		  for (int i =0; i< H->nb_sons();++i){
		    Block<double>& Hk = H->get_son(i);
		    if ((Hk.get_target_cluster().get_offset()<= oft) and (Hk.get_target_cluster().get_size())>=szt){
			for (int j =0; j< K->nb_sons();++j){
			  Block<double>& Kk= K->get_son(j);
			  if  ((Kk.get_source_cluster().get_offset()<= ofs) and (Kk.get_source_cluster().get_size())>=szs){
			    if((Kk.get_target_cluster().get_offset() == Hk.get_source_cluster().get_offset()) and (Kk.get_target_cluster().get_size()==Hk.get_source_cluster().get_size())){
			      vtemp.push_back(&Hk); vtemp.push_back(&Kk); }	} }}}}}
		      //cout<< "bb"<< endl;
            //normalement vtemp c'est notre nouveaux Sh mais il a peut être des low rank
            //On parcours et on regarde si il en a a des low rank;

            int nn = vtemp.size()/2;
		cout<<"!!!!!!"<< vtemp.size()<< ','<< nn << endl;
            //cout<<"!!!!!!!!!!!!!!!!!!"<<nn<<endl;
            for (int l = 0; l< nn ; ++l){
cout<< "aie"<< endl;
                Block<double>* H = vtemp[2*l]; 
cout<< "aieaie"<< endl;
		Block<double>* K= vtemp[2*l+1];
cout<< "aieaieaie" << endl;
            // on regarde si les deux sont low rank
                if(!(H->get_low_rank_block_data() ==nullptr) and !(K->get_low_rank_block_data() == nullptr)){
                    Matrix<double> Uh,Vh,Uk,Vk;
                    Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
                    Uk = K->get_low_rank_block_data()->Get_U();Vk = K->get_low_rank_block_data()->Get_V();
                    Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
		    cout<<"c'est quoi le pb:"<< V.nb_rows()<<','<<V.nb_cols()<<endl;
                    int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
                    sr.push_back(U); sr.push_back(V); of.push_back(rt); of.push_back(rs);
		cout<<"ici"<< endl;
                }
                //Celle de droite low rank
                else if( !(K->get_low_rank_block_data() == nullptr)){
		    cout<<"la"<< endl;
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
	            cout<< "ou" << endl;
                    Matrix<double> U = H->get_low_rank_block_data()->Get_U();
                    Matrix<double> V = H->get_low_rank_block_data()->Get_V();
                    Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());

                    for (int rr = 0 ; rr <V.nb_rows();++rr){
cout<< V.nb_rows() <<','<< rr<<':'<<"ok" << endl;
                        const vector<double> x = V.get_row(rr);

                        vector<double> y (K->get_target_cluster().get_size(),0.);
cout<< "ok" << endl;
                        produitP(K,x,y,'T');
                        for(int kl = 0 ; kl< y.size(); ++kl){
                            W(rr,kl)=y[kl];}
                    }
cout<<"héhé"<< endl;
                    sr.push_back(U);sr.push_back(W); of.push_back(H->get_target_cluster().get_offset()); of.push_back(K->get_source_cluster().get_offset());
cout<< "ok ?" << endl;
                }
                else{
		cout<< "lala"<< endl;
                    //ni H ni K ne sont low rank on les laisse dans vh
                    vh.push_back(H);vh.push_back(K);}
cout<<"l,nn="<<l<<','<<nn<< endl;
            }


            SumExpression Res(sr,of,vh);
            cout<<"on sort de retrict"<< endl;
            return Res;}
	
void tri(){
	vector<Matrix<double>> sr = this->get_sr();
	vector<Block<double>*> vtemp = this->get_sh();
	vector<Block<double>*> vh; vector<int> of = this->get_off();
	int N  =  vtemp.size()/2; //cout<< N<< endl;
	for (int k =0; k < N; ++k){
		//cout<< k << endl;
		Block<double>* H = vtemp[2*k]; Block<double>* K = vtemp[2*k+1];
		if( (H->get_low_rank_block_data() !=nullptr) and (K->get_low_rank_block_data() != nullptr)){
		    //cout<<"ici"<< endl;
                    Matrix<double> Uh,Vh,Uk,Vk;
                    Uh = H->get_low_rank_block_data()->Get_U();Vh = H->get_low_rank_block_data()->Get_V();
                    Uk = K->get_low_rank_block_data()->Get_U();Vk = K->get_low_rank_block_data()->Get_V();
                    Matrix<double> U = Uh; Matrix<double> V = Vh*(Uk*Vk);
		    //cout<<"c'est quoi le pb:"<< V.nb_rows()<<','<<V.nb_cols()<<endl;
                    int rt = H->get_target_cluster().get_offset(); int rs = K->get_source_cluster().get_offset();
                    sr.push_back(U); sr.push_back(V); of.push_back(rt); of.push_back(rs);
                }
                //Celle de droite low rank
                else if( !(K->get_low_rank_block_data() == nullptr)){
		    //cout<<"la"<< endl;
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
	            //cout<< "ou" << endl;
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
                else {
		//cout<< "lala"<< endl;
                    //ni H ni K ne sont low rank on les laisse dans vh
                    vh.push_back(H);vh.push_back(K);}
            }
//cout<< "ohohohohhohohhooé" << endl;
	this->SR= sr; this->SH = vh; this->off = of;
		}
/*
   //fonction qui va être appelée sur les feuilles non admissibles
     Matrix<double> Evaluate(const int& m,const int&  n){
	cout<< "on est rentré"<< endl;
     vector<Matrix<double>> Sr= this->get_sr(); vector<Block<double>*> Sh = this->get_sh();
     Matrix<double> res(m,n);
     for (int k =0; k < (Sr.size())/2; ++k){
       Matrix<double> U = Sr[2*k]; Matrix<double> V = Sr[2*k+1];
	//cout<<"+++++++++++++++++++"<< U(0,0) << endl;
       res = res +U*V;}
	cout<<"il a fait les low rank" << endl; 
     for (int k =0; k< (Sh.size())/2; ++k){
	//cout<< k << endl;
       Block<double>* H = Sh[2*k]; Block<double>* K = Sh[2*k+1];
	cout<<"#######################"<<endl;
	cout<< H->nb_sons() <<endl;
	cout<< K->nb_sons()<< endl;

	if(H->get_dense_block_data()== nullptr){
	cout<< "00000"<< endl;}
	if(H->get_low_rank_block_data() == nullptr){
	cout<< "11111"<< endl; }
	if(K->get_dense_block_data()== nullptr){
	cout<< "222222"<< endl;}
	if(K->get_low_rank_block_data() == nullptr){
	cout<< "3"<< endl; }

	cout<< " pour H " << H->get_target_cluster().get_nb_sons() << ','<< H->get_source_cluster().get_nb_sons()<<endl;
	cout<< "pour K "  << K->get_target_cluster().get_nb_sons() << ','<< K->get_source_cluster().get_nb_sons()<< endl;
	cout<<"------------"<< H->get_target_cluster().get_son(0).get_size() << endl;

       DenseBlockData<double> hh = *H->get_dense_block_data();
       DenseBlockData<double> kk = *K->get_dense_block_data();
       Matrix<double> hk = (hh)* (kk);
       res= res+hk;}
	cout << "wtf" << endl;
     return res;}

*/
     Matrix<double> Evaluate(const int& m,const int&  n){
	cout<<"Evaluate" << endl;
     vector<Matrix<double>> Sr= this->get_sr(); vector<Block<double>*> Sh = this->get_sh();
     Matrix<double> res(m,n);
     // On évalue les rangs faibles
     for (int k =0; k < (Sr.size())/2; ++k){
       Matrix<double> U = Sr[2*k]; Matrix<double> V = Sr[2*k+1];
       res = res +U*V;}
	cout<< "rang faibles ok" << endl;
     //On évalue les blocs
     // Attention apparament il se peut que SH ne contiennent pas que des feuilles du coup faut remettre un coup de réucursion et évaluer les sous blocs ( si il y en a ) et puis les additionner au bloc cortrespondant de la matrice résultat
     for (int k =0; k< (Sh.size())/2; ++k){
       Block<double>* H = Sh[2*k]; Block<double>* K = Sh[2*k+1];
       // cas facile ou H et K ne sont pas des Hmat
       if (H->nb_sons() == 0 and K->nb_sons() ==0 ){
       	  DenseBlockData<double> hh = *H->get_dense_block_data();
          DenseBlockData<double> kk = *K->get_dense_block_data();
          Matrix<double> hk = (hh)* (kk);
          res= res+hk;}
      // cas où il faut continuer à descendre, Il va faloir differencier tout les cas ...
      else{
		// cas hmat* full -> hmat*vecteur
		if ( (H->nb_sons() > 0) and (K->nb_sons() ==0) ){
	   	   Matrix<double> U = *K->get_dense_block_data();
			cout<< U.nb_cols()<< endl;
                    Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
                    for (int rr = 0 ; rr <U.nb_cols();++rr){
                        const vector<double> x = U.get_col(rr);
                        int rp = H->get_target_cluster().get_size();
                        vector<double> y (rr,0.);
                        produitP(H,x,y,'N');
                        for(int kl = 0 ; kl< y.size(); ++kl){
                            W(kl,rr)=y[kl];}}
}
	  } }
     return res;}
   

 };
//classe pour avoir un virtualgenerator a partir d'une matrice
class MyMatrix2 : public VirtualGenerator<double> {
  int space_dim;
  Matrix<double> mat;

  public:
    // Constructor
  MyMatrix2(int space_dim0, int nr, int nc, const Matrix<double> M) : VirtualGenerator(nr, nc), space_dim(space_dim0),mat(M) {}

    // Virtual function to overload
    double get_coef(const int &k, const int &j) const {
      return mat(k,j);
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
 //pettit getter


    int get_rows(){ return p1.size();}
    int get_cols(){ return p2.size();}

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


/*

void Hmult(Block<double>*L,Block<double>* Lts, SumExpression HK, const VirtualCluster& t,const VirtualCluster& s,const double& epsilon, const int& spacedim = 2){
  // comment on fait si tsons=0 et pas s? on va pas rentrer dans la boucle =========> on est obligé de différencier les cas ?? ====> et la réponse est oui evidemment ........ C'est completement con ca veut dire qu'on doit faire les trois cas 
  //recursion: si on est pas une feuillle on découpe et on recommence sur le fils avec la bonne restriction
  int tt = t.get_nb_sons(); int ss = s.get_nb_sons();
  if( (tt>0) and (ss>0)) {
	cout<<"yo"<< endl;
    int rep =0;
    for(int k = 0; k< t.get_nb_sons(); ++k){
      for(int l = 0; l< s.get_nb_sons(); ++l){
	Lts->build_son(t.get_son(k),s.get_son(l));
	Block<double>& Lrep = Lts->get_son(rep);
	SumExpression restr = HK.Restrict(t.get_son(k),s.get_son(l));
	Hmult(L,&Lrep,restr,t.get_son(k),s.get_son(l),epsilon,spacedim);
	rep+=1;}
    } }
  else if( (ss > 0) ){
	cout<< "yoyo"<< endl;
	for( int k =0 ; k < ss ; ++k){
		Lts->build_son(t,s.get_son(k));
		Block<double>& Lrep = Lts->get_son(k);
		SumExpression restr = HK.Restrict(t,s.get_son(k));
		Hmult(L,&Lrep,restr,t,s.get_son(k),epsilon,spacedim);} }
  else if( tt>0){
	cout<< "yoyoyo" << endl;
	for( int k =0 ; k < tt ; ++k){
		Lts->build_son(t.get_son(k),s);
		Block<double>& Lrep = Lts->get_son(k);
		SumExpression restr = HK.Restrict(t.get_son(k),s);
		Hmult(L,&Lrep,restr,t.get_son(k),s,epsilon,spacedim);} }
		
  else{
    cout<< "aa"<< ss<<','<< tt<< endl;
    if (Lts->IsAdmissible()){
	cout<<"bb"<<endl;
      // ici il faudrait faire la troncature mais pour l'instant on va juste sommer
      // Normalement il fauddrait faire res= lowrank(res+UV)
      vector<Matrix<double>> sr = HK.get_sr();
      int r = sr[0].nb_cols();
      Matrix<double> val (Lts->get_target_cluster().get_size(),Lts->get_source_cluster().get_size());
	//LowRankMatrix < double> val( sr[0],sr[1]);
      for( int i = 0 ; i < sr.size()/2 ; ++ i ){
	//LowRankMatrix<double> temp (sr[2*i],sr[2*i+1]);
	val = val+ sr[2*i]*sr[2*i+1];
//val= val +temp;
}

      const Block<double>& bb = *Lts;
      MyMatrix2 gen (spacedim ,val.nb_rows(),val.nb_cols(),val);
      //const VirtualLowRankGenerator<double> &LRGenerator = sympartialACA<double>;
      //Lts-> LowRankGenerator = make_shared<sympartialACA<double>>;
	cout<<"cc"<<endl;
      auto LRgenerator= *new sympartialACA<double>;
      //Lts->get_low_rank_block_data() 
      auto ll= new LowRankMatrix<double>(*Lts, gen, LRgenerator, 0, 0, r, epsilon, true);
      //new LowRankMatrix<double> lr(*Lts, gen, make_shared<sympartialACA<double>>(), 0, 0, r, epsilon, true);
      //Lts->get_block_data() = std::unique_ptr<VirtualBlockData<double>>(Lts->get_low_rank_block_data());
      auto ls = Lts->get_block_data();
      unique_ptr<VirtualBlockData<double>> temp (ll);
      auto tt = Lts->get_block_data();
      tt= &*temp;  // est ce que ca suffit pour modfier L puisque en soit Lts et un bloc de L?
      //temp->swap(ll);
      //ls->swap(unique_ptr(ll);
      //ls = unique_ptr<VirtualBlockData<double>>(ll);
      //Block<double>* brb = *L->Get_Tasks();
      //brb.push_back(Lts);
//vector<Block<double>*>* BB= &brb;
  //    &L->get_tasks() = *BB;
}
    else{
cout<<"ou alors"<< endl;
	cout<< HK.get_sr().size()<<','<< HK.get_sh().size()<< endl;
	HK.tri();
	cout<< HK.get_sr().size()<<','<< HK.get_sh().size()<< endl;
      Matrix<double> val = HK.Evaluate(Lts->get_target_cluster().get_size(), Lts->get_source_cluster().get_size() );
      //cout<< val(0,0)<<','<< val(1,1)<< ',' << val(2,2) << endl;
	//cout<<1<< endl;
	cout<<"ici"<<endl;
	MyMatrix2 gen ( spacedim,val.nb_rows(),val.nb_cols(), val);
	//cout<< 2<< endl;
      //MyMatrix2 gen(spacedim, val.nb_rows(),val.nb_cols() );
	//cout<< val.nb_rows() <<','<< val.nb_cols()<< endl;
	//cout<< gen.norm() << endl;
	cout<< "la"<< endl;
	auto ll = new DenseBlockData<double>(*Lts, gen,true);
	//cout<< 3<< endl;
	cout<<"ici la"<< endl;	
	auto ls  = Lts->get_block_data();
	//cout<< 4<< endl;
cout<<"mais pas"<< endl;
	unique_ptr<VirtualBlockData<double>> temp(ll);
	//cout<< 5 << endl;
cout<<"?"<< endl;
	ls = &*temp;
cout<<"!"<< endl;
      //Lts->get_dense_block_data() = new DenseBlockData<double>(*Lts, &val, true);
      //Lts->block_data       = std::unique_ptr<VirtualBlockData<double>>(L->dense_block_data);
    //L->get_tasks()->push_back(Lts);
}
  }
}
*/
void Hmult2(Block<double>*L,Block<double>* Lts, SumExpression HK, const VirtualCluster& t,const VirtualCluster& s,const double& epsilon, const int& spacedim = 2){
  // comment on fait si tsons=0 et pas s? on va pas rentrer dans la boucle =========> on est obligé de différencier les cas ?? ====> et la réponse est oui evidemment ........ C'est completement con ca veut dire qu'on doit faire les trois cas 
  //recursion: si on est pas une feuillle on découpe et on recommence sur le fils avec la bonne restriction
  int tt = t.get_nb_sons(); int ss = s.get_nb_sons();cout<< "t,s" << ' '<< tt<<','<<ss<< endl;
  if( (tt>0) or (ss>0)) {
	cout<<"yo"<< endl;
    if ( (tt>0) and (ss>0) ){
    int rep =0;
    for(int k = 0; k< t.get_nb_sons(); ++k){
      for(int l = 0; l< s.get_nb_sons(); ++l){

	Lts->build_son(t.get_son(k),s.get_son(l));
	Block<double>& Lrep = Lts->get_son(rep);
	cout<< t.get_son(k).get_size() <<',' << s.get_son(l).get_size() << endl;  
	SumExpression restr = HK.Restrict(t.get_son(k),s.get_son(l));
	cout<< "test" << endl;
	Hmult2(L,&Lrep,restr,t.get_son(k),s.get_son(l),epsilon,spacedim);
	rep+=1;}
    } ;} 
   if( (ss > 0) and (tt==0)){
	cout<< "yoyo"<< endl;
	for( int k =0 ; k < ss ; ++k){
		Lts->build_son(t,s.get_son(k));
		Block<double>& Lrep = Lts->get_son(k);
		SumExpression restr = HK.Restrict(t,s.get_son(k));
		Hmult2(L,&Lrep,restr,t,s.get_son(k),epsilon,spacedim);} }
  else if( (tt>0) and (ss ==0) ) {
	cout<< "yoyoyo" << endl;
	for( int k =0 ; k < tt ; ++k){
		Lts->build_son(t.get_son(k),s);
		Block<double>& Lrep = Lts->get_son(k);
		SumExpression restr = HK.Restrict(t.get_son(k),s);
		Hmult2(L,&Lrep,restr,t.get_son(k),s,epsilon,spacedim);} }
		}
  else{
    cout<< "aa"<< ss<<','<< tt<< endl;
    if (Lts->IsAdmissible()){
	cout<<"bb"<<endl;
      // ici il faudrait faire la troncature mais pour l'instant on va juste sommer
      // Normalement il fauddrait faire res= lowrank(res+UV)
      vector<Matrix<double>> sr = HK.get_sr();
      int r = sr[0].nb_cols();
      Matrix<double> val (Lts->get_target_cluster().get_size(),Lts->get_source_cluster().get_size());
      for( int i = 0 ; i < sr.size()/2 ; ++ i ){
	val = val+ sr[2*i]*sr[2*i+1];
}
      const Block<double>& bb = *Lts;
      MyMatrix2 gen (spacedim ,val.nb_rows(),val.nb_cols(),val);
      //const VirtualLowRankGenerator<double> &LRGenerator = sympartialACA<double>;
      //Lts-> LowRankGenerator = make_shared<sympartialACA<double>>;
	cout<<"cc"<<endl;
      auto LRgenerator= *new sympartialACA<double>;
      //Lts->get_low_rank_block_data() 
      auto ll= new LowRankMatrix<double>(*Lts, gen, LRgenerator, 0, 0, r, epsilon, true);
      //new LowRankMatrix<double> lr(*Lts, gen, make_shared<sympartialACA<double>>(), 0, 0, r, epsilon, true);
      //Lts->get_block_data() = std::unique_ptr<VirtualBlockData<double>>(Lts->get_low_rank_block_data());
      auto ls = Lts->get_block_data();
      unique_ptr<VirtualBlockData<double>> temp (ll);
      auto tt = Lts->get_block_data();
      tt= &*temp;  // est ce que ca suffit pour modfier L puisque en soit Lts et un bloc de L?
      //temp->swap(ll);
      //ls->swap(unique_ptr(ll);
      //ls = unique_ptr<VirtualBlockData<double>>(ll);
      //Block<double>* brb = *L->Get_Tasks();
      //brb.push_back(Lts);
//vector<Block<double>*>* BB= &brb;
  //    &L->get_tasks() = *BB;
}
    else{
// Le problème c'est qu'on appelle evaluate alors que dans la sum expresion il y a encore des H et K bloc et du coup les deux ont un BlockData = nullptr ce qui fait planter Evaluate

// Plusieurs solutions: - Changer la descente au début de Hmult pour être sure / changer Evaluate pour qu'il continue à descendre si on lyui donne deux blocs
cout<<"ou alors"<< endl;
	if (Mycondition().ComputeAdmissibility(t,s,1)) {
 		cout << 454545454545454<< endl; }
	cout<< t.get_nb_sons()<<','<< s.get_nb_sons()<< endl;
	cout<< HK.get_sr().size()<<','<< HK.get_sh().size()<< endl;
	HK.tri();
	cout<< HK.get_sr().size()<<','<< HK.get_sh().size()<< endl;
      Matrix<double> val = HK.Evaluate(Lts->get_target_cluster().get_size(), Lts->get_source_cluster().get_size() );
	cout<<"ici"<<endl;
	MyMatrix2 gen ( spacedim,val.nb_rows(),val.nb_cols(), val);
	cout<< "la"<< endl;
	auto ll = new DenseBlockData<double>(*Lts, gen,true);
	//cout<< 3<< endl;
	cout<<"ici la"<< endl;	
	auto ls  = Lts->get_block_data();
cout<<"mais pas"<< endl;
	unique_ptr<VirtualBlockData<double>> temp(ll);
cout<<"?"<< endl;
	ls = &*temp;
cout<<"!"<< endl;
}
  } }



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
    double eta     = 1;

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

     cout<< "Test Restrict"<< endl;
     pair<vector<Matrix<double>>,vector<int>> sr; vector<Block<double>*> sh;
     Block<double>& B= Rt->get_son(0);Block<double>* BB = &B;
     sh.push_back(BB); sh.push_back(BB);
     pair<pair<vector<Matrix<double>>,vector<int>>,vector<Block<double>*>> S;
     S.first = sr; S.second=sh;
     SumExpression AB(BB,BB);
     SumExpression Stest =  AB.Restrict(BB->get_target_cluster().get_son(0),BB->get_source_cluster().get_son(1));
     vector<Matrix<double>> srtest = Stest.get_sr(); vector<Block<double>*> ssh = Stest.get_sh(); vector<int> offs = Stest.get_off();
     cout<< srtest.size()<<','<<offs.size()<<endl;
     cout<< ssh.size()<<endl;
     cout<< ssh[0]->get_target_cluster().get_offset()<<','<<ssh[0]->get_target_cluster().get_size()<<endl;
     cout<<'!'<< BB->get_target_cluster().get_son(0).get_offset()<<','<<BB->get_source_cluster().get_son(1).get_offset();
     cout<<'!'<< BB->get_target_cluster().get_son(0).get_size()<<','<<BB->get_source_cluster().get_son(1).get_size()<<endl;

     for ( int k =0; k< 2; ++k){
         cout<<"++++++++++++++++"<<endl;
         cout<< ssh[2*k]->get_target_cluster().get_offset()<<','<<ssh[2*k]->get_target_cluster().get_size()<<endl;
         cout<< ssh[2*k]->get_source_cluster().get_offset()<<','<<ssh[2*k]->get_source_cluster().get_size()<<endl;
         cout<<"------------------"<< endl;
         cout<< ssh[2*k+1]->get_target_cluster().get_offset()<<','<<ssh[2*k+1]->get_target_cluster().get_size()<<endl;
         cout<< ssh[2*k+1]->get_source_cluster().get_offset()<<','<<ssh[2*k+1]->get_source_cluster().get_size()<<endl;

     }
     // --------------------> OK ? faut tester avec evaluate, le problème c'est que evaluate on peut l'utiliser que sur des feuilels
// liste de low rank
vector<LowRankMatrix<double>> lowrank ;
vector<DenseBlockData<double>> full ;
vector<Block<double>*> Full;
vector<Block<double>*> Lr;
for (int k =0; k< Data.size(); ++k){
	Block<double>* br = Data[k];
	if( br->get_low_rank_block_data() == nullptr ) {
		auto temp = br->get_dense_block_data();
		full.push_back(*temp ) ;
		Full.push_back(br);}
	else{
		auto  temp = br->get_low_rank_block_data();
		lowrank.push_back( *temp );
		Lr.push_back(br);}
}

cout<< full.size() << endl;
cout<< lowrank.size() << endl;

LowRankMatrix<double> l0 = lowrank[0];
Matrix<double> uu =  l0.Get_U();
cout<< uu(0,0) <<','<< uu(1,0) << uu(2,0)<< endl;
DenseBlockData<double> d0 = full[0];
cout << d0.nb_rows()<< ','<< d0.nb_cols()<< endl;
cout<< d0(0,0)<< endl;
//////
Block<double>* ff= Full[0];
SumExpression s (ff,ff);

Matrix<double> test = s.Evaluate(37,37);
cout<< test(0,0) << endl;
Matrix<double> ppp = d0*d0;
cout<<"!!!!!!!!!!!!!!!!!!!!"<< endl;
cout<< normFrob(test-ppp) << endl;
cout<< normFrob(test)<<";;"<< normFrob(ppp) << endl;
// Evaluate marche bien sur les blocs full
vector<int> ofs; 
//ofs.push_back(0); ofs.push_back(0);
vector<Block<double>*> sm; 
sm.push_back(ff);sm.push_back(ff);sm.push_back(ff);sm.push_back(ff);
vector<Matrix<double>> ss ; 
SumExpression ttest(ss,ofs,sm);

Matrix<double> tpt = ttest.Evaluate(37,37);
Matrix<double> ref = d0*d0 + d0*d0;
cout<< "/*/*/*/*/*"<<normFrob(tpt - ref)<< endl;
// Ca à l'aire de marcher

// Test pour descendre dans le block
Block<double> L(&*make_shared<Mycondition>(), *HA.get_target_cluster(),*HA.get_source_cluster() );
L.set_eta(1);
/*
Split(L);
int rep = 0;
Block<double>* pt;
pt = &L;
while(pt->nb_sons() >  0){
	rep+=1;
	pt= &pt->get_son(0);}
cout<< rep<< endl;
*/
// -----------> Ca descend bien dans la hiérarchie en créant les fils
// On va essayer traquer Restrict sur qq itérations
// On initialise notre sumexpr = A*A -> sr=[], sh=[bloc A,bloc A]
SumExpression se(&B,&B);
/*
cout<< B.get_target_cluster().get_nb_sons()<< endl;
for (int k =0; k< B.get_target_cluster().get_nb_sons(); ++k){
	for(int l=0; l < B.get_source_cluster().get_nb_sons(); ++l){
		SumExpression stemp = se.Restrict( B.get_target_cluster().get_son(k), B.get_source_cluster().get_son(l) );
		cout<<"1)"<< stemp.get_sr().size()<<','<<se.get_sh().size()<< endl;
		cout<<"2)"<< stemp.get_sh()[0]->get_target_cluster().get_offset()<< endl;
		cout<<"3)"<< stemp.get_sh()[0]->get_source_cluster().get_offset()<< endl;
		cout<< "4)"<<stemp.get_sh()[1]->get_target_cluster().get_offset()<< endl;
		cout<< "5)"<<stemp.get_sh()[1]->get_source_cluster().get_offset()<< endl;} }

bool test0 = true;
int aa = 0;
while(aa<8){
	cout<<"boulcle :"<<aa<< endl;
	se = se.Restrict(se.get_sh()[0]->get_target_cluster().get_son(0),se.get_sh()[1]->get_source_cluster().get_son(0));
cout<< "nb son"<< se.get_sh()[0]->get_target_cluster().get_nb_sons() << endl;
   	//test = (se.get_sh()[0]->get_target_cluster().get_son(0).get_nb_sons() > 0);
	cout<< "#####################r"<< endl;
	for(int k =0 ; k< se.get_sh().size()/2; ++k){
	cout<< se.get_sh()[2*k]->get_target_cluster().get_offset()<<','<< se.get_sh()[2*k]->get_source_cluster().get_offset()<<"_____________"<<se.get_sh()[2*k]->get_target_cluster().get_size()<<','<< se.get_sh()[2*k]->get_source_cluster().get_size()<< endl;
	cout<< se.get_sh()[2*k+1]->get_target_cluster().get_offset()<<','<< se.get_sh()[2*k+1]-> get_source_cluster().get_offset() <<"_____________"<<se.get_sh()[2*k+1]->get_target_cluster().get_size()<<','<< se.get_sh()[2*k+1]->get_source_cluster().get_size()<< endl;}

	aa+=1;
}
cout<< se.get_sr().size()<< ',' << se.get_sh().size()<< endl;
Matrix<double> bloc = se.Evaluate(se.get_sh()[0]->get_target_cluster().get_size(), se.get_sh()[1]->get_source_cluster().get_size() ) ;
Matrix<double> m(HA.get_target_cluster()->get_size(),HA.get_source_cluster()->get_size());
     for (int k =0 ;k < HA.get_target_cluster()->get_size();++k){
       for (int l =0; l< HA.get_source_cluster()->get_size();++l){
	 m(k,l) = A.get_coef(k,l);
       }}
     Matrix<double> mm = m*m;
 
Matrix<double> vra ( se.get_sh()[0]->get_target_cluster().get_size(), se.get_sh()[1]->get_source_cluster().get_size());
for (int k = 0 ; k < se.get_sh()[0]->get_target_cluster().get_size() ; ++k){
	for ( int l = 0; l< se.get_sh()[1]->get_source_cluster().get_size() ; ++l){
		vra(k,l) = mm(k,l) ; } } 

cout<<"tadaaaaaaaa:"<< normFrob(bloc-vra)/normFrob(vra) << endl;
*/
/*
for(int k =0 ; k< se.get_sh().size()/2; ++k){
	cout<< se.get_sh()[2*k]->get_target_cluster().get_offset()<<','<< se.get_sh()[2*k]->get_source_cluster().get_offset()<<"_____________"<<se.get_sh()[2*k]->get_target_cluster().get_size()<<','<< se.get_sh()[2*k]->get_source_cluster().get_size()<< endl;
	cout<< se.get_sh()[2*k+1]->get_target_cluster().get_offset()<<','<< se.get_sh()[2*k+1]-> get_source_cluster().get_offset() <<"_____________"<<se.get_sh()[2*k+1]->get_target_cluster().get_size()<<','<< se.get_sh()[2*k+1]->get_source_cluster().get_size()<< endl;}
*/


/*
     cout<<"test evaluate"<< endl;

     // On va regarder dans A des blocs non admissibles carré et on va essyer de regarder si ca marche
     // C'est sensé nous renvoyé le produit de deux matrices pleines
     vector<Block<double>*> test;
     for (int k =0 ; k < Data.size(); ++k){
       Block<double>* b = Data[k];
       if (b->get_target_cluster().get_size() == b->get_source_cluster().get_size()){
	 if ( !(b->get_dense_block_data() == nullptr) ){
	   test.push_back(b);} } }
     cout<<test.size()<< endl;

     Block<double>* b = test[0];
     int ot = b->get_target_cluster().get_offset(); int os = b->get_source_cluster().get_offset();
     int st = b->get_target_cluster().get_size() ; int ss = b->get_source_cluster().get_size();
     Matrix<double> aa (st,ss);
     for (int i = 0 ; i < st; ++i){
	 for(int j =0; j < ss; ++j){
	   aa(i,j)= A.get_coef(i+ot,j+os);}
     }
     DenseBlockData<double> bb = *b->get_dense_block_data();
     Matrix<double> ma = bb*bb;
     SumExpression seval(b,b);
     Matrix<double> mtest = seval.Evaluate(st,ss);
     cout<< normFrob(mtest-ma)/normFrob(ma)<< endl;
     // for (int r=0; r<15; ++r){
     //   cout<<mtest(1,r)<<','<<ma(1,r)<<endl;}
     // Ca marche pas quand on test contre A mais est ce que c'est à cause des permutations
*/


     //test hmult//
     cout<< "test Hmult"<< endl;
/*
     Block<double> L(Mycondition, HA.get_target_cluster(),HA.get_source_cluster() );
     Matrix<double> m(HA.get_target_cluster().get_size(),HA.get_source_cluster().get_size());
     for (int k =0 ;k < HA.get_target_cluster().get_size()){
       for (int l =0; l< HA.get_source_cluster().get_size()){
	 m(k,l) = A.get_coef(k,l);
       }}
     Matrix<double> mm = m*m;
     SumExpression AA(Rt,Rt);
     Hmult(L,AA,HA.get_target_cluster(),HA.get_source_cluster(),epsilon, 2);
     cout << Frobenius_absolute_error(L, mm) / normFrob(mm) <<endl;
    */   
	   
   
	//Block<double> L(&*make_shared<Mycondition>(), *HA.get_target_cluster(),*HA.get_source_cluster() );
	Matrix<double> m(HA.get_target_cluster()->get_size(),HA.get_source_cluster()->get_size());
        for (int k =0 ;k < HA.get_target_cluster()->get_size(); ++k){
            for (int l =0; l< HA.get_source_cluster()->get_size(); ++l){
	        m(k,l) = A.get_coef(k,l);
            }}
        Matrix<double> mm = m*m;
        SumExpression AA(Rt,Rt);
        Hmult2(&L,&L,AA,*HA.get_target_cluster(),*HA.get_source_cluster(),epsilon, 2);
        //cout << Frobenius_absolute_error(L, mm) / normFrob(mm) <<endl;



     

MPI_Finalize();
}
