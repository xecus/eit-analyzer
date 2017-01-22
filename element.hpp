class Element{
	private:
	public:
	Eigen::MatrixXd pos;
	Eigen::VectorXd fai;
	Eigen::VectorXd psi;	
	int point_num[8];
	int NumOfPoint;
	int NumOfDim;
	~Element(){
	}
	Element(){
	}

	//Init
	void InitDim2(){
		this->NumOfPoint = 4;
		this->NumOfDim = 2;
		pos = Eigen::MatrixXd::Zero(4,2);
		fai = Eigen::VectorXd::Zero(4);
		psi = Eigen::VectorXd::Zero(4);
	}
	void InitDim3(){
		this->NumOfPoint = 8;
		this->NumOfDim = 3;
		pos = Eigen::MatrixXd::Zero(8,3);
		fai = Eigen::VectorXd::Zero(8);
		psi = Eigen::VectorXd::Zero(8);
	}

	//J Matrix
	Eigen::Matrix<double,2,2>  Dim2_JMatrix(double n1,double n2){
		Eigen::Matrix<double,2,2> res;
		res = Eigen::MatrixXd::Zero(2,2);
		for(int i=0;i<4;i++){
			res(0,0) += Dim2_Nx(i,0,n1,n2) * pos(i,0);
			res(0,1) += Dim2_Nx(i,0,n1,n2) * pos(i,1);
			res(1,0) += Dim2_Nx(i,1,n1,n2) * pos(i,0);
			res(1,1) += Dim2_Nx(i,1,n1,n2) * pos(i,1);
		}
		return res;
	}
	Eigen::Matrix<double,3,3> Dim3_JMatrix(double n1,double n2,double n3){
		Eigen::Matrix<double,3,3> res;
		res = Eigen::MatrixXd::Zero(3,3);
		for(int i=0;i<8;i++){
			res(0,0) += Dim3_Nx(i,0,n1,n2,n3) * pos(i,0);
			res(0,1) += Dim3_Nx(i,0,n1,n2,n3) * pos(i,1);
			res(0,2) += Dim3_Nx(i,0,n1,n2,n3) * pos(i,2);
			res(1,0) += Dim3_Nx(i,1,n1,n2,n3) * pos(i,0);
			res(1,1) += Dim3_Nx(i,1,n1,n2,n3) * pos(i,1);
			res(1,2) += Dim3_Nx(i,1,n1,n2,n3) * pos(i,2);
			res(2,0) += Dim3_Nx(i,2,n1,n2,n3) * pos(i,0);
			res(2,1) += Dim3_Nx(i,2,n1,n2,n3) * pos(i,1);
			res(2,2) += Dim3_Nx(i,2,n1,n2,n3) * pos(i,2);
		}
		return res;
	}

	//det(J Matrix)
	double Dim2_JMatrixVal(double n1,double n2){
		Eigen::Matrix<double,2,2> J = Dim2_JMatrix(n1,n2);
		return J(0,0)*J(1,1) - J(1,0) * J(0,1);
	}
	double Dim3_JMatrixVal(double n1,double n2,double n3){
		Eigen::Matrix<double,3,3> J = Dim3_JMatrix(n1,n2,n3);
		double JMatrix_determinant = 0;
		JMatrix_determinant += J(0,0)*J(1,1)*J(2,2);
		JMatrix_determinant += J(0,1)*J(1,2)*J(2,0);
		JMatrix_determinant += J(0,2)*J(1,0)*J(2,1);
		JMatrix_determinant -= J(0,2)*J(1,1)*J(2,0);
		JMatrix_determinant -= J(0,1)*J(1,0)*J(2,2);
		JMatrix_determinant -= J(0,0)*J(1,2)*J(2,1);
		return JMatrix_determinant;
	}

	//Nx2
	Eigen::Matrix<double,2,4> Dim2_Nx2(double n1,double n2){
		Eigen::Matrix<double,2,2> J = Dim2_JMatrix(n1,n2);
		Eigen::FullPivLU< Eigen::Matrix<double,2,2> > lu( J );
		Eigen::Matrix<double,2,2> J_Inverse =  lu.inverse();
		Eigen::Matrix<double,2,4> temp;
		temp = Eigen::MatrixXd::Zero(2,4);
		for(int i=0;i<4;i++){
			temp(0,i) = Dim2_Nx(i,0,n1,n2);
			temp(1,i) = Dim2_Nx(i,1,n1,n2);
		}
		return J_Inverse * temp;
	}
	Eigen::Matrix<double,3,8> Dim3_Nx2(double n1,double n2,double n3){
		Eigen::Matrix<double,3,3> J = Dim3_JMatrix(n1,n2,n3);
		Eigen::FullPivLU< Eigen::Matrix<double,3,3> > lu( J );
		Eigen::Matrix<double,3,3> J_Inverse =  lu.inverse();
		Eigen::Matrix<double,3,8> temp;
		temp = Eigen::MatrixXd::Zero(3,8);
		for(int i=0;i<8;i++){
			temp(0,i) += Dim3_Nx(i,0,n1,n2,n3);
			temp(1,i) += Dim3_Nx(i,1,n1,n2,n3);
			temp(2,i) += Dim3_Nx(i,2,n1,n2,n3);
		}
		return J_Inverse * temp;
	}

	//K Matrix(Origin)
 	Eigen::Matrix<double,4,4> Dim2_KMatrix_origin(double n1,double n2){
 		Eigen::Matrix<double,2,4> temp = Dim2_Nx2(n1,n2);
 		Eigen::Matrix<double,4,4> prod = temp.transpose() * temp;
 		return prod;
	}
	Eigen::Matrix<double,8,8> Dim3_KMatrix_origin(double n1,double n2,double n3){
		Eigen::Matrix<double,3,8> temp = Dim3_Nx2(n1,n2,n3);
		Eigen::Matrix<double,8,8> prod = temp.transpose() * temp;
		return prod;
	}

	//K Matrix
 	Eigen::Matrix<double,4,4> Dim2_KMatrix(){
 		Eigen::Matrix<double,4,4> res;
 		res = Eigen::MatrixXd::Zero(4,4);
		for(int i=0;i<gi.bunten_num;i++){
 			for(int j=0;j<gi.bunten_num;j++){
 				double temp = gi.weight[i] * gi.weight[j]
				* Dim2_JMatrixVal(gi.bunten[i],gi.bunten[j]);
 				res += Dim2_KMatrix_origin(gi.bunten[i],gi.bunten[j]) * temp;
 			}
 		}
 		return res;
	}
 	Eigen::Matrix<double,8,8> Dim3_KMatrix(){
		Eigen::Matrix<double,8,8> res;
		res = Eigen::MatrixXd::Zero(8,8);
		for(int i=0;i<gi.bunten_num;i++){
			for(int j=0;j<gi.bunten_num;j++){
				for(int k=0;k<gi.bunten_num;k++){
					double temp = gi.weight[i] * gi.weight[j] * gi.weight[k]
					* Dim3_JMatrixVal(gi.bunten[i],gi.bunten[j],gi.bunten[k]);
					res += Dim3_KMatrix_origin(gi.bunten[i],gi.bunten[j],gi.bunten[k]) * temp;
				}
			}
		}
		return res;
	}

	//CalcFaiDimm
	double Dim2_CalcFaiDiff_fai(int nn,double n1,double n2){
		double res=0.0;
		for(int i=0;i<4;i++)
			res+=Dim2_Nx(i,nn,n1,n2) * fai(i);
		return res;
	}
	double Dim2_CalcFaiDiff_psi(int nn,double n1,double n2){
		double res=0.0;
		for(int i=0;i<4;i++)
			res+=Dim2_Nx(i,nn,n1,n2) * psi(i);
		return res;
	}
	double Dim3_CalcFaiDiff_fai(int nn,double n1,double n2,double n3){
		double res=0.0;
		for(int i=0;i<8;i++)
			res+=Dim3_Nx(i,nn,n1,n2,n3) * fai(i);
		return res;
	}
	double Dim3_CalcFaiDiff_psi(int nn,double n1,double n2,double n3){
		double res=0.0;
		for(int i=0;i<8;i++)
			res+=Dim3_Nx(i,nn,n1,n2,n3) * psi(i);
		return res;
	}

	//Nx
	double Dim2_Nx(int n,double n1,double n2){
		double res = 1 / 4.0;
		if(n==0) res*=(1-n1)*(1-n2);
		if(n==1) res*=(1+n1)*(1-n2);
		if(n==2) res*=(1+n1)*(1+n2);
		if(n==3) res*=(1-n1)*(1+n2);
		return res;
	}
	double Dim3_Nx(int n,double n1,double n2,double n3){
		double res;
		res =1.0/8.0;
		if(n==0) res *= (1-n1)*(1-n2)*(1-n3);
		if(n==1) res *= (1+n1)*(1-n2)*(1-n3);
		if(n==2) res *= (1+n1)*(1+n2)*(1-n3);
		if(n==3) res *= (1-n1)*(1+n2)*(1-n3);
		if(n==4) res *= (1-n1)*(1-n2)*(1+n3);
		if(n==5) res *= (1+n1)*(1-n2)*(1+n3);
		if(n==6) res *= (1+n1)*(1+n2)*(1+n3);
		if(n==7) res *= (1-n1)*(1+n2)*(1+n3);
		return res;
	}

	//Nx
	double Dim2_Nx(int n,int nn,double n1,double n2){
		double res = 1 / 4.0;
		if(n==0){
			if(nn==0) res*=(-1)*(1-n2);
			if(nn==1) res*=(1-n1)*(-1);
		}
		if(n==1){
			if(nn==0) res*=(1)*(1-n2);
			if(nn==1) res*=(1+n1)*(-1);
		}
		if(n==2){
			if(nn==0) res*=(1)*(1+n2);
			if(nn==1) res*=(1+n1)*(1);
		}
		if(n==3){
			if(nn==0) res*=(-1)*(1+n2);
			if(nn==1) res*=(1-n1)*(1);
		}
		return res;
	}
	double Dim3_Nx(int n,int nn,double n1,double n2,double n3){
		double res;
		res =1.0/8.0;
		if(n==0){
			if(nn==0) res *= (-1)*(1-n2)*(1-n3);
			if(nn==1) res *= (1-n1)*(-1)*(1-n3);
			if(nn==2) res *= (1-n1)*(1-n2)*(-1);
		}
		if(n==1){
			if(nn==0) res *= (1)*(1-n2)*(1-n3);
			if(nn==1) res *= (1+n1)*(-1)*(1-n3);
			if(nn==2) res *= (1+n1)*(1-n2)*(-1);
		}
		if(n==2){
			if(nn==0) res *= (1)*(1+n2)*(1-n3);
			if(nn==1) res *= (1+n1)*(1)*(1-n3);
			if(nn==2) res *= (1+n1)*(1+n2)*(-1);
		}
		if(n==3){
			if(nn==0) res *= (-1)*(1+n2)*(1-n3);
			if(nn==1) res *= (1-n1)*(1)*(1-n3);
			if(nn==2) res *= (1-n1)*(1+n2)*(-1);
		}
		if(n==4){
			if(nn==0) res *= (-1)*(1-n2)*(1+n3);
			if(nn==1) res *= (1-n1)*(-1)*(1+n3);
			if(nn==2) res *= (1-n1)*(1-n2)*(1);
		}
		if(n==5){
			if(nn==0) res *= (1)*(1-n2)*(1+n3);
			if(nn==1) res *= (1+n1)*(-1)*(1+n3);
			if(nn==2) res *= (1+n1)*(1-n2)*(1);
		}
		if(n==6){
			if(nn==0) res *= (1)*(1+n2)*(1+n3);
			if(nn==1) res *= (1+n1)*(1)*(1+n3);
			if(nn==2) res *= (1+n1)*(1+n2)*(1);
		}
		if(n==7){
			if(nn==0) res *= (-1)*(1+n2)*(1+n3);
			if(nn==1) res *= (1-n1)*(1)*(1+n3);
			if(nn==2) res *= (1-n1)*(1+n2)*(1);
		}
		return res;
	}

};
