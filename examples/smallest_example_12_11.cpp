#include <htool/htool.hpp>
#include <memory>
#include <iomanip>
using namespace std;
using namespace htool;


// Condition d'adm 
class Mycondition final : public VirtualAdmissibilityCondition {
  public:
    bool ComputeAdmissibility(const VirtualCluster &target, const VirtualCluster &source, double eta) const override {
        bool admissible = 1 * std::min(target.get_rad(), source.get_rad()) < eta * std::max((norm2(target.get_ctr() - source.get_ctr()) - target.get_rad() - source.get_rad()), 0.);
        return admissible; }
};

// Produit hmat/bloc vecteur

void produitP(Block<double>& B, const vector<double>& x, vector<double>& y, char op){
    int of_t = B.get_target_cluster().get_offset(); int of_s = B.get_source_cluster().get_offset();
    int sz_t = B.get_target_cluster().get_size();  int sz_s = B.get_source_cluster().get_size();
    if (!(B.get_block_data() ==nullptr) and (B.nb_sons()==0)){
        if (op =='N'){ B.get_block_data()->add_mvprod_row_major(x.data()+of_s,y.data()+of_t,1,'N',op); }
        else{ B.get_block_data()->add_mvprod_row_major(x.data()+of_t,y.data()+of_s,1,'N',op);} }
    else{
        for (int r =0; r<B.nb_sons();++r){
            Block<double>& Br = B.get_son(r); produitP(Br,x,y,op);} }
}

//classe pour avoir un virtualgenerator a partir d'une matrice
class MyMatrix2 : public VirtualGenerator<double> {
  int space_dim;
  Matrix<double> mat;

  public:
    // Constructor
  MyMatrix2(int space_dim0, int nr, int nc, const Matrix<double> M) : VirtualGenerator(nr, nc), space_dim(space_dim0),mat(M) {}

    // Virtual function to overload
    double get_coef(const int &k, const int &j) const {
      return mat(k,j);}

    // Virtual function to overload
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const override {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                ptr[j + M * k] = this->get_coef(rows[j], cols[k]);}}}

    // Matrix vector product
    std::vector<double> operator*(std::vector<double> a) {
        std::vector<double> result(nr, 0);
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                result[j] += this->get_coef(j, k) * a[k]; }}
        return result;}

    // Frobenius norm
    double norm() {
        double norm = 0;
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                norm += this->get_coef(j, k); } }
        return norm; }
};

//classe Sum-expr
 class SumExpression {
  private :
   // liste de matrice de rang faible U1,V1,U2,V2..., liste de leurs offset, liste de hmat H1,K1,H2,K2
    vector<Matrix<double>> SR;
    vector<int> off;
    vector<Block<double>> SH;

  public:
  
    //Constructeur qui prend deux hmat : L=HK = sumexpr(I,I)
    SumExpression( Block<double>& H,  Block<double>& K){ vector<Block<double>> vb; vb.push_back(H); vb.push_back(K);}
  /*
// Constructeur avec tout
    SumExpression(const vector<Matrix<double>>& SR0,const vector<int>& off0,const vector<Block<double>*>& Sh0){this->SR=SR0;this->off=off0;this->SH=Sh0;}

   //getters
    vector<Matrix<double>> get_sr(){return this->SR;}
    vector<int> get_off(){return this->off;}
    vector<Block<double>*> get_sh(){return this->SH;}

    //Fonction Tri -> vérifie qu'il y a pas de low rank qui trainent dans sh
    // Il faudrait mieux le faire direct dans restrict puisque ce rajoute un O(#S) 
    //-> une fonction qui étant donné H,K (au moins 1 lr) le plug dans sr et off
    void Tri(){
	vector<Block<double>*> sh = this->get_sh(); vector<Matrix<double>> sr= this->get_sr(); 
	vector<Block<double>*> shres;
	vector<int> of = this->get_off();int nh = sh.size()/2; 
	for (int i =0 ; i< nh ; ++i){
		Block<double>* H = sh[2*i]; Block<double>* K = sh[2*i+1];
		int oft = H->get_target_cluster().get_offset(); int ofs = K->get_source_cluster().get_offset();
		
		//On différencie lr*lr, lr*hmat, hmat*lr
		//lrmat lrmat
		if( (H->get_low_rank_block_data() != nullptr) and (K->get_low_rank_block_data() != nullptr) ) {
			Matrix<double> Uh = H->get_low_rank_block_data()->Get_U(); Matrix<double> Vh = H->get_low_rank_block_data()->Get_V();
			Matrix<double> Uk = K->get_low_rank_block_data()->Get_U(); Matrix<double> Vk = K->get_low_rank_block_data()->Get_V();
			Matrix<double> U = Uh* ( Vh * Uk); sr.push_back(U); sr.push_back(Vk);
		    of.push_back(oft); of.push_back(ofs); }
		
		// hmat lrmat -> plein de hmat vecteus
		if ( ( H->get_low_rank_block_data()== nullptr) and (K->get_low_rank_block_data() != nullptr) ) {
			Matrix<double> U = K->get_low_rank_block_data()->Get_U(); Matrix<double> V = K->get_low_rank_block_data()->Get_V();
                        Matrix<double> W( H->get_target_cluster().get_size(),U.nb_cols());
                        for (int j = 0 ; j <U.nb_cols();++j){
                        	const vector<double> x = U.get_col(j); vector<double> y (H->get_target_cluster().get_size(),0.); 
				produitP(*H,x,y,'N');
                        	for(int kl = 0 ; kl< y.size(); ++kl){ W(kl,j)=y[kl];}}
                        sr.push_back(W);sr.push_back(V); of.push_back(oft); of.push_back(ofs);}
		
		// lrmat hmat -> (hmat*vect)^T
		if ( (H->get_low_rank_block_data() != nullptr) and ( K->get_low_rank_block_data() == nullptr) ){
			Matrix<double> U = H->get_low_rank_block_data()->Get_U(); Matrix<double> V = H->get_low_rank_block_data()->Get_V();
			Matrix<double> W( U.nb_cols(),K->get_source_cluster().get_size());
                        for (int rr = 0 ; rr <V.nb_rows();++rr){
                        	const vector<double> x = V.get_row(rr); vector<double> y (K->get_target_cluster().get_size(),0.);
                        	produitP(*K,x,y,'T');
                        	for(int kl = 0 ; kl< y.size(); ++kl){ W(rr,kl)=y[kl];} }
                    	sr.push_back(U);sr.push_back(W); of.push_back(oft); of.push_back(ofs); }
		
		// les deux sont hmat on les laisse dans sh
		if ( (H->get_low_rank_block_data()== nullptr) and (K->get_low_rank_block_data() == nullptr) ){
			shres.push_back(H); shres.push_back(K); } 
		}
		SR = sr;SH = shres; off = of; }
		

    //Fonction Restrict
    SumExpression Restrict(const VirtualCluster& t, const VirtualCluster& s){
		vector<Matrix<double>> sr =  this->SR; vector<Block<double>*> sh = this->SH; vector<int> of = this->off;
		vector<Matrix<double>> sres; vector<Block<double>*> shres; vector<int> ofres;
		int oft = t.get_offset();  int szt = t.get_size(); int ofs = s.get_offset() ; int szs = s.get_size(); 
		int nr = sr.size()/2 ; int nh = sh.size()/2; 
	
	//On commence par prendre les restrictions des low rank
	for (int n = 0 ; n < nr ; ++n){
		Matrix<double> U= sr[2*n]; Matrix<double> V = sr[2*n+1]; int ofu = of[2*n] ; int ofv = of[2*n+1];
		
		// On doit restreinde les lignes de U à t et les colonnes de V à s
		// Cette opération n'est possible que si tau est dans les indices des lignes de U et s dans les indices des colonnnes de V
		if ( (ofu >= oft) and ( szt+(oft-ofu) <= U.nb_rows() )  and ( ofv >= ofs ) and ( szs+(ofs-ofv))){
			Matrix<double> Uk( szt,U.nb_cols() ); Matrix<double> Vk ( V.nb_rows(), szs);
			for (int i =0 ; i < szt ; ++i){
				for (int j =0 ; j< U.nb_cols() ; ++j){ Uk(i,j) = U(oft-ofu + i, j); } }
			for (int i=0; i < V.nb_rows() ; ++i ){
				for (int j =0; j < szs; ++j){ Vk(i,j)=V(i,ofs-ofv+j); } }
			sres.push_back(Uk); sres.push_back(Vk); ofres.push_back(oft); ofres.push_back(ofs) ; } }
	
	//On restreint maintenant les hmat grace à leurs fils
	// on part du princpe qu'on ne peut pas restreindre sur des éléments qui ne sont pas dans l'arbre
	for (int n =0; n < nh; ++n){
		//cout<<"!!"<<n <<endl;
		Block<double>* H=  sh[2*n]; Block<double>* K = sh[2*n+1];
		
		// On regarde si par hasard ils ont pas déja la même taille auquel cas on les met direct dans shres
		if ( (H->get_target_cluster().get_size() == szt) and ( H->get_target_cluster().get_offset() == oft) and ( K->get_source_cluster().get_size() ==szs) and (K->get_source_cluster().get_offset() == ofs ) ){
			shres.push_back(H); shres.push_back(K); } 
		
		//Sinon on prend les bon fils de H et K
		else {
			for (int k =0; k< H->nb_sons(); ++k){
				Block<double>& Hk = H->get_son(k) ; 
				// pas besoin de tester si H a pas le bon offset et taille ( on suppose qu'on peut restreindre qu'au fils)
				int szh = Hk.get_target_cluster().get_size(); int ofh = Hk.get_target_cluster().get_offset();
				int szr = Hk.get_source_cluster().get_size(); int ofr = Hk.get_source_cluster().get_offset();
				if ((szt== szh) and ( ofh==oft)){
					for (int l =0; l< K->nb_sons() ; ++l){
						Block<double>& Kk = K->get_son(l);
						int szk = Kk.get_source_cluster().get_size(); int ofk= Kk.get_source_cluster().get_offset();
						int szrr = Kk.get_target_cluster().get_size(); int ofrr= Kk.get_target_cluster().get_offset();
						// On regarde si K a le bon source et si H.source = K.target
						if( (szk == szs) and (ofk == ofs) and (szr== szrr) and (ofr ==ofrr) ){
							shres.push_back(&Hk); shres.push_back(&Kk);
							//--------->  C'est ici qu'il faudrait faire le tri
							} } } } }
	
}
	//Il se peut que dans notre shres il yait des low rank don il faut la trier
	SumExpression res(sres, ofres, shres); res.Tri();
	return res;}
	
	//Fonction evaluate
	
	Matrix<double> Evaluate(const int& m, const int & n){
		Matrix<double> Res(m,n); vector<Block<double>*> sh= this->get_sh(); vector<Matrix<double>> sr = this->get_sr();
		//les low rank
		for (int k =0; k < sr.size()/2; ++k){
			Matrix<double> U= sr[2*k]; Matrix<double> V = sr[2*k+1];
			Res = Res+ U*V;}
		//Les blocs ( normalement quand on l'appelle c'est sur des feuilles non admissibles et donc de taille min
		for (int k =0; k < sh.size()/2 ; ++k ){
			Block<double>* H = sh[2*k]; Block<double>* K = sh[2*k+1];
			Matrix<double> temp = *H->get_dense_block_data() * *K->get_dense_block_data() ;
			Res = Res +temp;}
		return Res; }
		
	// Fonction troncature
	
	LowRankMatrix<double> Troncature(Block<double>* B,int rank, double epsilon, bool use_perm){
		this->Tri(); 
		Matrix <double> t0 = this->get_sr()[0]* this->get_sr()[1] + this->get_sr()[2]* this->get_sr()[3];
		MyMatrix2 T0(3, t0.nb_rows(), t0.nb_cols(), t0);
		LowRankMatrix<double> T (*B,T0, *make_shared<sympartialACA<double>>(),0,0,rank, epsilon,use_perm);
		for (int k =2; k< this->get_sr().size()/2; ++k){
			Matrix<double> tk = T.Get_U()*T.Get_V()+ this->get_sr()[2*k] * this->get_sr()[2*k+1] ; 
			MyMatrix2 Tk(3, tk.nb_rows(), tk.nb_cols(), tk);
			LowRankMatrix<double> t (*B,Tk, *make_shared<sympartialACA<double>>() ,0,0,rank, epsilon,use_perm);
			T = t;	}
		return T;}
*/		
	

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
        return (1.) / (4 * M_PI * std::sqrt(1e-5 + std::inner_product(p1.begin() + space_dim * k, p1.begin() + space_dim * k + space_dim, p2.begin() + space_dim * j, double(0), std::plus<double>(), [](double u, double v) { return (u - v) * (u - v); }))); }

    // Virtual function to overload
    void copy_submatrix(int M, int N, const int *const rows, const int *const cols, double *ptr) const override {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                ptr[j + M * k] = this->get_coef(rows[j], cols[k]);} } }
 //pettit getter
    int get_rows(){ return p1.size();}
    int get_cols(){ return p2.size();}
    // Matrix vector product
    std::vector<double> operator*(std::vector<double> a) {
        std::vector<double> result(nr, 0);
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
                result[j] += this->get_coef(j, k) * a[k];}}
        return result;}

    // Frobenius norm
    double norm() {
        double norm = 0;
        for (int j = 0; j < nr; j++) {
            for (int k = 0; k < nc; k++) {
      	        norm += this->get_coef(j, k); }}
        return norm;}
};

//Fonction Hmult
// on veut calculer L=H*K du coup faudrait un void dans Block.hpp qui agisse sur un bloc vide instancier avant et le bloc courant ce serait this
// Pour l'instant on va juste faire un Proto


/*

void Hmult(Block<double>* L, Block<double>* Lk , SumExpression S,int rank, double epsilon, bool use_perm) {
	// Si Lt n'est pas un feuille (comment on vérifie BIEN ca?) on construit ses fils et on rappel Hmult avec la sumexpr restreint au fils
	//Récursion
	//cout<< Lk->get_target_cluster().get_offset()<<','<< Lk->get_target_cluster().get_size()<< endl;
	//cout << Lk -> get_source_cluster().get_offset()<< ',' << Lk->get_source_cluster().get_size()<< endl;
	cout << "+++++++++++" << endl; 
	cout<< Lk->get_target_cluster().get_nb_sons() << ',' << Lk->get_source_cluster().get_nb_sons()<< endl;
	cout<< " ++++++++" << endl; 
	
	if( !((Lk->get_target_cluster().get_nb_sons() == 0) and (Lk->get_source_cluster().get_nb_sons() == 0)) ){
		//On est obligé de différencier les trois cas ...
		if ( (Lk->get_target_cluster().get_nb_sons() > 0 ) and ( Lk->get_source_cluster().get_nb_sons() > 0) ){
			cout<< "1"<< endl;
			int repere = 0; 
			for( int k =0 ; k< Lk->get_target_cluster().get_nb_sons() ; ++k ){ 
				for (int l =0; l< L->get_source_cluster().get_nb_sons(); ++l){
					Lk->build_son(Lk->get_target_cluster().get_son(k),Lk->get_source_cluster().get_son(l));
					SumExpression Skl = S.Restrict(Lk->get_target_cluster().get_son(k),Lk->get_source_cluster().get_son(l));
					Hmult(L,&Lk->get_son(repere),Skl, rank, epsilon , use_perm);  repere +=1;}}
		cout<< "1 ok" << endl;}
		else if ( (Lk->get_target_cluster().get_nb_sons() > 0 )){
			cout<< "2" << endl; 
			int repere =0;
		for( int k =0 ; k< Lk->get_target_cluster().get_nb_sons() ; ++k ){ 
				Lk->build_son(Lk->get_target_cluster().get_son(k),Lk->get_source_cluster());
				SumExpression Skl = S.Restrict(Lk->get_target_cluster().get_son(k),Lk->get_source_cluster());
				Hmult(L,&Lk->get_son(repere),Skl,rank,epsilon , use_perm);  repere +=1;}
				cout<< "2 ok" << endl;} 
		else{ cout<< "3"<< endl;
			int repere =0;
			for (int l =0; l< L->get_source_cluster().get_nb_sons(); ++l){
					Lk->build_son(Lk->get_target_cluster(),Lk->get_source_cluster().get_son(l));
					SumExpression Skl = S.Restrict(Lk->get_target_cluster(),Lk->get_source_cluster().get_son(l));
					Hmult(L,&Lk->get_son(repere),Skl,rank, epsilon, use_perm);  repere +=1;}
					cout<< "3 ok" << endl;}
		}
	else{cout<< "yo"<< endl;
		//Normalement on est sur un feuille du coup on va devoir passer le *block_data en !nullptr 
		//Premier cas feuille non admissible -> on fait Evaluate sur la sumexpr et on fait pointer le dense block data dessus
		if ( !(Lk->IsAdmissible()) )  { 
		    cout<< "a" << endl;
			Matrix<double> Val = S.Evaluate(Lk->get_target_cluster().get_size(),Lk->get_source_cluster().get_size());
			cout<< "a0" << endl;
			MyMatrix2 val (3, Val.nb_rows(),Val.nb_cols(), Val);
			cout<< Val(0,0) << ','<< Val(1,0) << endl;
			cout<< "a1" << endl;
			if ( Lk->get_block_data() == nullptr){cout<< "yoyoyo" << endl;}
			Lk->compute_dense_block(val, use_perm);
			cout<< "a2" << endl; 
}
		// Deuxieme cas on est admissible du coup normalement on a que Sh=[] ( on va quand même faire un tri)
		// On fait pointer le bloc data sur la troncature
		else { 
			cout<< "b" << endl;
			LowRankMatrix<double> T = S.Troncature(Lk, rank , epsilon , use_perm) ; 
			MyMatrix2 gen (3, T.nb_rows(), T.nb_cols(), T.Get_U() * T.Get_V() );
			Lk->compute_low_rank_block(rank,epsilon,gen,*make_shared<sympartialACA<double>>() ,0,0, use_perm );

			}
		}
		cout<< "?!"<< endl;
	}

*/

// Il nous faut une fonction pour calculer la norm de frob d'un bloc
/*
double norm ( Block<double>* L, Block<double>* Lk , Matrix<double > A){
	if ( Lk->get_nb_sons() == 0 ) {
		auto val = Lk->get_block_data(); int oft = Lk->get_target_cluster().get_offset();
		int szt = Lk->get_target_cluster()get_size(); 
		int ofs = Lk->get_source_target_cluster().get_offset() ; int szs = Lk->get_source_cluster().get_size();
		for  } }
		*/
		
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
    double epsilon = 0.001;	double eta     = 1;

    // n² points on a regular grid in a square
    int n    = std::sqrt(4761);	int size = n * n;

    // p1: points in a square in the plane z=z1
    double z = 1;
    vector<double> p(3 * size);
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            p[3 * (j + k * n) + 0] = j;	p[3 * (j + k * n) + 1] = k;	p[3 * (j + k * n) + 2] = z;	}
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

    HMatrix<double> Res(t, t, epsilon, eta, 'S', 'U');
    Res.build(A, p.data());
    vector<Block<double>*> Data = HA.get_ComputedBlock();vector<Block<double>*> DataRes = Res.get_ComputedBlock();
    //On va faire un liste de low rank et de full
    vector<Block<double>*> LowR, Full ;
    for (int k =0; k< DataRes.size() ; ++k ){
    	Block<double>* b = DataRes[k];
    	if( b->get_dense_block_data()== nullptr){ 		LowR.push_back(b); if(b->get_low_rank_block_data()==nullptr){cout << "wtf"<< endl;} }
    	else {		Full.push_back(b);}
    	}
	cout<< LowR.size()<< ','<< Full.size()<<','<<LowR.size()+Full.size()-DataRes.size()<< endl;
    Block<double>* bb0 = Data[0]; Block<double>* Rt = bb0->get_root();
    // Je sais pas pourquoi mais la root à qu'un fils, elle même, ca casse un peu tout du coup je part de son(0)
    Block<double>* L = &Rt->get_son(0);
    
    /////////////////////
     cout <<"Test Produit"<<endl;
     vector<double> xx(4761,2), yy(4761,0), ref(4761,0),xx_perm(4761,0),yy_perm(4761,0);
     global_to_cluster(&Rt->get_source_cluster(),xx.data(),xx_perm.data());
     cluster_to_global(&Rt->get_target_cluster(),yy.data(),yy_perm.data());
     produitP(*Rt,xx_perm,yy_perm,'T');
     cluster_to_global(&Rt->get_target_cluster(),yy_perm.data(),yy.data());
     vector<double> yref = A*xx;
     cout<< norm2(yy-yref)/norm2(yref)<<endl;
     cout<< endl;
     //--------------------> Ok
     //////////////////////
    vector<Block<double>> vv;
    
    /* 
     
     /////////////////////////
    cout<< " test sumexpr"<< endl;
    
    SumExpression S(L,L);
    cout<< S.get_sr().size() << ',' << S.get_off().size()<< ',' << S.get_sh().size()<<endl;
    //-----------------
    cout<<"test Tri:"<< endl<<endl;
    vector<Block<double>*> v; v.push_back(L);v.push_back(L);v.push_back(LowR[1]);v.push_back(LowR[2]);
    vector<Matrix<double>> lr; vector<int> ii;
    SumExpression Stri(lr,ii,v);
    cout<< Stri.get_sr().size()<< ','<< Stri.get_off().size()<<','<< Stri.get_sh().size()<< endl;
    Stri.Tri();
    cout<< Stri.get_sr().size()<< ','<< Stri.get_off().size()<<','<< Stri.get_sh().size()<< endl;
    cout<< endl<<endl;
    //------------------> Ok
	cout<< "Test Restrict"<< endl;
	// on restreint (L,L) 
	cout<<"on retreint à t(0),s(0)="<< L->get_son(0).get_target_cluster().get_offset()<<','<<L->get_son(0).get_source_cluster().get_offset() << endl;
	cout<< "de taille " << L->get_son(0).get_target_cluster().get_size()<<','<<L->get_son(0).get_source_cluster().get_size()<< endl;
	SumExpression Srestr = S.Restrict(L->get_son(0).get_target_cluster(),L->get_son(0).get_source_cluster());
	cout<< Srestr.get_sr().size()<<','<<Srestr.get_off().size()<<','<<Srestr.get_sh().size()<< endl;
	cout<< "Les blocs de ce produits correspondent aux cluter: oft,ofs*oft,ofs   szt,szs"<<endl;
	cout<<'('<<Srestr.get_sh()[0]->get_target_cluster().get_offset()<<','<<Srestr.get_sh()[0]->get_source_cluster().get_offset()<<')'<<'x'<<'('<< Srestr.get_sh()[1]->get_target_cluster().get_offset()<<','<<Srestr.get_sh()[1]->get_source_cluster().get_offset()<<')'<< endl;
	cout<<'('<<Srestr.get_sh()[0]->get_target_cluster().get_size()<<','<<Srestr.get_sh()[0]->get_source_cluster().get_size()<<')'<<'x'<<'('<< Srestr.get_sh()[1]->get_target_cluster().get_size()<<','<<Srestr.get_sh()[1]->get_source_cluster().get_size()<<')'<< endl;
	cout<<"++++++++++++++++++++++"<<endl;
	cout<<'('<<Srestr.get_sh()[2]->get_target_cluster().get_offset()<<','<<Srestr.get_sh()[2]->get_source_cluster().get_offset()<<')'<<'x'<<'('<< Srestr.get_sh()[3]->get_target_cluster().get_offset()<<','<<Srestr.get_sh()[3]->get_source_cluster().get_offset()<<')'<< endl;
	cout<<'('<<Srestr.get_sh()[2]->get_target_cluster().get_size()<<','<<Srestr.get_sh()[2]->get_source_cluster().get_size()<<')'<<'x'<<'('<< Srestr.get_sh()[3]->get_target_cluster().get_size()<<','<<Srestr.get_sh()[3]->get_source_cluster().get_size()<<')'<< endl;
	
	//------------------> offset et taille ok
	cout<< "test evaluate"<<endl;
	Block<double>* bb =  Full[0]; cout<< bb->get_target_cluster().get_size() << ',' << bb->get_source_cluster().get_size()<< endl;
	vector<Block<double>*> vb; vb.push_back(bb); vb.push_back(bb); vb.push_back(bb); vb.push_back(bb); 
	SumExpression btest (lr,ii,vb); 
	Matrix<double> M = btest.Evaluate(74,74);
	cout<< normFrob( M - *bb->get_dense_block_data() * *bb->get_dense_block_data()- *bb->get_dense_block_data() * *bb->get_dense_block_data())/normFrob(*bb->get_dense_block_data() * *bb->get_dense_block_data() + *bb->get_dense_block_data() * *bb->get_dense_block_data() ) << endl;
	
	
	//-----------------> Ok
	int rmax = 100;
	cout<< "test Hmult" << endl;
	Block<double> HK ( &*make_shared<Mycondition> (), L->get_target_cluster(), L->get_source_cluster() );
	Hmult(&HK,&HK,S, rmax, epsilon , true);
	//Matrix<double> AA = A*A;
	//cout<< Frobenius_absolute_error(HK, AA) << endl ; 
	
*/
	
	

//////



    


     

MPI_Finalize();
}
