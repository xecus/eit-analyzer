
class CRSConverter{
private:
public:

	int NumOfDimM;
	int NumOfDimN;
	int NumOfMat;

	std::vector<int>		I_dat;
	std::vector<int>		J_dat;
	std::vector<double>		val_lhs_dat;
	std::vector<double>		val_rhs_dat;

	int 					*tmp_I_dat;
	int 					*tmp_J_dat;
	double					*tmp_val_lhs_dat;

	int 					*out_I_dat;
	int 					*out_J_dat;
	int 					out_nz;
	double 					*out_val_lhs_dat;
	double 					*out_val_rhs_dat;

	~CRSConverter(){
		Clear();
	}
	CRSConverter(){
		tmp_I_dat = NULL;
		tmp_J_dat = NULL;
		tmp_val_lhs_dat = NULL;
		Clear();
		return;
	}

	void Synth(){
		out_I_dat = &(I_dat.begin()[0]);
		out_J_dat = &(I_dat.begin()[0]);
		out_nz = J_dat.size();
		out_val_lhs_dat = &(val_lhs_dat.begin()[0]);
		out_val_rhs_dat = &(val_rhs_dat.begin()[0]);
		return;
	}

	void Clear(){
		if(tmp_I_dat != NULL)
			delete[] tmp_I_dat;
		if(tmp_J_dat != NULL)
			delete[] tmp_J_dat;
		if(tmp_val_lhs_dat != NULL)
			delete[] tmp_val_lhs_dat;
		I_dat.clear();
		J_dat.clear();
		val_lhs_dat.clear();
		val_rhs_dat.clear();
		NumOfDimM = 0;
		NumOfDimN = 0;
		NumOfMat = 0;
		return;
	}

	void Show(){
		std::cout << "I_dat=" << I_dat.size() << std::endl;
		std::cout << "J_dat=" << J_dat.size() << std::endl;
		std::cout << "val_lhs_dat=" << val_lhs_dat.size() << std::endl;
		std::cout << "val_rhs_dat=" << val_rhs_dat.size() << std::endl;
		std::cout << "NumOfMat=" << NumOfMat << std::endl;
		int sum=0;
		sum += I_dat.size() * sizeof(int);
		sum += J_dat.size() * sizeof(int);
		sum += val_lhs_dat.size() * sizeof(double);
		sum += val_rhs_dat.size() * sizeof(double);
		std::cout << "Size:" << sum << "Byte" << std::endl;
		std::cout << "Size:" << sum/(1024) << "KByte" << std::endl;
		std::cout << "Size:" << sum/(1024*1024) << "MByte" << std::endl;
		std::cout << "Size:" << sum/(1024*1024*1024) << "GByte" << std::endl;
		/*
		int *p_int;
		double *p_double;
		p_int = &(I_dat.begin()[0]);
		p_int = &(J_dat.begin()[0]);
		p_double = &(val_lhs_dat.begin()[0]);
		p_double = &(val_rhs_dat.begin()[0]);
		*/
		return;
	}

	void Add(
		Eigen::SparseMatrix<double,Eigen::RowMajor> A
		,Eigen::SparseMatrix<double,Eigen::RowMajor> B
		){

		int M,N,nz;
		M = A.rows();
		N = A.cols();
		nz = A.nonZeros();
		//std::cout << "(Add) M=" << M << " N=" << N << " nz=" << nz << std::endl;
		if( M != N ){
			std::cout << "(Error) Can't Add Matrix(M!=N)" << std::endl;
			return;
		}
		if( NumOfMat == 0 ){	//Init
			NumOfDimM = M;
			NumOfDimN = N;
			if (tmp_I_dat==NULL)	tmp_I_dat	=	new int[ N + 1 ];	//row-data(compressed)
			if (tmp_J_dat==NULL)	tmp_J_dat	=	new int[ nz ];		//col-data
			if (tmp_val_lhs_dat==NULL)	tmp_val_lhs_dat	=	new double[ nz ];	//val-data
		}else{
			if(NumOfDimM != M || NumOfDimN != N){
				std::cout << "(Error) Can't Add Matrix(Dim)" << std::endl;
				return;
			}
		}
		//non-zero element
		int cnt1=0,cnt2=0,prow = -1;
		for(int i=0;i<A.outerSize();++i){
			for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(A,i);
			it;++it){
				int row = it.row();
				int col = it.col();
				double val = it.value();
				tmp_val_lhs_dat[cnt1] = (double)val;
				tmp_J_dat[cnt1] = col;
				if( prow != row ){
					tmp_I_dat[cnt2++] = cnt1;
					prow = row;
				}
				cnt1++;
			}
		}
		tmp_I_dat[cnt2++] = cnt1;

		/*
		std::cout << "[I (" << (N+1) << ") ]" << std::endl;
		for(int i=0;i< N+1; i++) std::cout << tmp_I_dat[i] << ",";
		std::cout << std::endl;
		std::cout << "[J (" << nz << ")]" << std::endl;
		for(int i=0;i< nz ; i++) std::cout << tmp_J_dat[i] << ",";
		std::cout << std::endl;
		std::cout << "[val (" << nz << ")]" << std::endl;
		for(int i=0;i< nz ; i++) std::cout << tmp_val_lhs_dat[i] << ",";
		std::cout << std::endl;
		*/

		//Generating(A)
		if(NumOfMat==0){
			for(int i=0;i< N+1; i++) I_dat.push_back( tmp_I_dat[i] );
		}else{
			int vback = I_dat.back();
			I_dat.pop_back();
			for(int i=0;i< N+1; i++) I_dat.push_back( tmp_I_dat[i] + vback );
		}
		for(int i=0;i< nz; i++)
			J_dat.push_back( tmp_J_dat[i] + NumOfMat * N );
		for(int i=0;i< nz; i++)
			val_lhs_dat.push_back( tmp_val_lhs_dat[i] );

		//Generating(B)
		for(int i=0;i<N;i++) val_rhs_dat.push_back(B.coeffRef(i,0));

		NumOfMat++;

		return;
	}


};