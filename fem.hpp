
class Fem{
private:
public:

	//Struct
	int SumOfElement;
	int SumOfPoint;
	Eigen::SparseMatrix<double> FEM_base;
	ESet *esp;		//Electrode Set
	Element *_ele;	//Element

	~Fem(){
		delete[] _ele;
	}

	Fem(ESet *_esp){

		this->esp = _esp;
		if( this->esp->adp->NumOfDim == 2 ) InitDim2();
		if( this->esp->adp->NumOfDim == 3 ) InitDim3();
		
		return;
	}

	void InitDim3(){
		this->SumOfElement = this->esp->adp->NumOfEleX * this->esp->adp->NumOfEleY * this->esp->adp->NumOfEleZ;
		this->SumOfPoint = (this->esp->adp->NumOfEleX+1)*(this->esp->adp->NumOfEleY+1)*(this->esp->adp->NumOfEleZ+1);
		double fem_h_x = this->esp->adp->LenghOfX / this->esp->adp->NumOfEleX;
		double fem_h_y = this->esp->adp->LenghOfY / this->esp->adp->NumOfEleY;
		double fem_h_z = this->esp->adp->LenghOfZ / this->esp->adp->NumOfEleZ;
		_ele = new Element[this->SumOfElement];
		for(int i=0;i<this->SumOfElement;i++) _ele[i].InitDim3();
		for(int z=0;z< this->esp->adp->NumOfEleZ;z++){
			for(int y=0;y< this->esp->adp->NumOfEleY;y++){
				for(int x=0;x< this->esp->adp->NumOfEleX;x++){
					int temp_num = this->esp->adp->ConvertElement(x,y,z);
					_ele[temp_num].pos(0,0)=fem_h_x*x;
					_ele[temp_num].pos(0,1)=fem_h_y*y;
					_ele[temp_num].pos(0,2)=fem_h_z*z;
					_ele[temp_num].pos(1,0)=fem_h_x*(x+1);
					_ele[temp_num].pos(1,1)=fem_h_y*y;
					_ele[temp_num].pos(1,2)=fem_h_z*z;
					_ele[temp_num].pos(2,0)=fem_h_x*(x+1);
					_ele[temp_num].pos(2,1)=fem_h_y*(y+1);
					_ele[temp_num].pos(2,2)=fem_h_z*z;
					_ele[temp_num].pos(3,0)=fem_h_x*x;
					_ele[temp_num].pos(3,1)=fem_h_y*(y+1);
					_ele[temp_num].pos(3,2)=fem_h_z*z;
					for(int i=0;i<4;i++){
						_ele[temp_num].pos(4+i,0)= _ele[temp_num].pos(i,0);
						_ele[temp_num].pos(4+i,1)= _ele[temp_num].pos(i,1);
						_ele[temp_num].pos(4+i,2)= _ele[temp_num].pos(i,2) + fem_h_z;
					}
					_ele[temp_num].point_num[0] = this->esp->adp->ConvertPoint(x,y,z);
					_ele[temp_num].point_num[1] = this->esp->adp->ConvertPoint(x+1,y,z);
					_ele[temp_num].point_num[2] = this->esp->adp->ConvertPoint(x+1,y+1,z);
					_ele[temp_num].point_num[3] = this->esp->adp->ConvertPoint(x,y+1,z);
					_ele[temp_num].point_num[4] = this->esp->adp->ConvertPoint(x,y,z+1);
					_ele[temp_num].point_num[5] = this->esp->adp->ConvertPoint(x+1,y,z+1);
					_ele[temp_num].point_num[6] = this->esp->adp->ConvertPoint(x+1,y+1,z+1);
					_ele[temp_num].point_num[7] = this->esp->adp->ConvertPoint(x,y+1,z+1);
				}
			}
		}
		
		return;
	}

	void InitDim2(){
		this->SumOfElement = this->esp->adp->NumOfEleX * this->esp->adp->NumOfEleY;
		this->SumOfPoint = (this->esp->adp->NumOfEleX+1)*(this->esp->adp->NumOfEleY+1);
		double fem_h_x = this->esp->adp->LenghOfX / this->esp->adp->NumOfEleX;
		double fem_h_y = this->esp->adp->LenghOfY / this->esp->adp->NumOfEleY;
		_ele = new Element[this->SumOfElement];
		for(int i=0;i<this->SumOfElement;i++) _ele[i].InitDim2();
		for(int y=0;y<this->esp->adp->NumOfEleY;y++){
			for(int x=0;x< this->esp->adp->NumOfEleX;x++){
				int temp_num = this->esp->adp->ConvertElement(x,y);
				_ele[temp_num].pos(0,0)=fem_h_x*x;
				_ele[temp_num].pos(0,1)=fem_h_y*y;
				_ele[temp_num].pos(1,0)=fem_h_x*(x+1);
				_ele[temp_num].pos(1,1)=fem_h_y*y;
				_ele[temp_num].pos(2,0)=fem_h_x*(x+1);
				_ele[temp_num].pos(2,1)=fem_h_y*(y+1);
				_ele[temp_num].pos(3,0)=fem_h_x*x;
				_ele[temp_num].pos(3,1)=fem_h_y*(y+1);
				_ele[temp_num].point_num[0] = this->esp->adp->ConvertPoint(x,y);
				_ele[temp_num].point_num[1] = this->esp->adp->ConvertPoint(x+1,y);
				_ele[temp_num].point_num[2] = this->esp->adp->ConvertPoint(x+1,y+1);
				_ele[temp_num].point_num[3] = this->esp->adp->ConvertPoint(x,y+1);
			}
		}
		return;
	}

	void SetSigma2D(Eigen::VectorXd &sigma){
		FEM_base.resize(SumOfPoint,SumOfPoint);
		for(int i=0;i<SumOfElement;i++){
			Eigen::Matrix<double,4,4> temp_k = _ele[i].Dim2_KMatrix();
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					double temp = temp_k(j,k) * sigma(i);
					FEM_base.coeffRef(_ele[i].point_num[j] , _ele[i].point_num[k]) += temp;
				}
			}
		}
		FEM_base.prune(0,0);
	}

	void SetSigma3D(Eigen::VectorXd &sigma){
		FEM_base.resize(SumOfPoint,SumOfPoint);
		for(int i=0;i<SumOfElement;i++){
			Eigen::Matrix<double,8,8> temp_k;
			//temp_k = _ele[i].Dim3_KMatrix();
			temp_k << 0.006128472222222218,0.001397569444444442,-0.0001345486111111116,0.001397569444444443
			,-0.002795138888888877,-0.002230902777777775,-0.00153211805555555,-0.002230902777777777
			,0.001397569444444442,0.006128472222222214,0.001397569444444443,-0.0001345486111111116
			,-0.002230902777777777,-0.002795138888888884,-0.002230902777777772,-0.001532118055555551
			,-0.0001345486111111116,0.001397569444444443,0.00612847222222221,0.001397569444444442
			,-0.001532118055555551,-0.002230902777777772,-0.002795138888888884,-0.002230902777777775
			,0.001397569444444443,-0.0001345486111111116,0.001397569444444442,0.00612847222222222
			,-0.002230902777777778,-0.001532118055555551,-0.002230902777777777,-0.002795138888888877
			,-0.002795138888888877,-0.002230902777777777,-0.001532118055555551,-0.002230902777777778
			,0.006128472222222221,0.001397569444444441,-0.0001345486111111119,0.001397569444444443
			,-0.002230902777777775,-0.002795138888888884,-0.002230902777777772,-0.001532118055555551
			,0.001397569444444441,0.006128472222222215,0.001397569444444443,-0.0001345486111111119
			,-0.00153211805555555,-0.002230902777777772,-0.002795138888888884,-0.002230902777777777
			,-0.0001345486111111119,0.001397569444444443,0.006128472222222213,0.001397569444444442
			,-0.002230902777777777,-0.001532118055555551,-0.002230902777777775,-0.002795138888888877
			,0.001397569444444443,-0.0001345486111111119,0.001397569444444442,0.006128472222222218;
			for(int j=0;j<8;j++){
				for(int k=0;k<8;k++){
					double temp = temp_k(j,k) * sigma(i);
					FEM_base.coeffRef(_ele[i].point_num[j] , _ele[i].point_num[k]) += temp;
				}
			}
		}
		FEM_base.prune(0,0);
	}

	/* Set Condition */
	void SetDirchletCondition(int point_num,double val
		,Eigen::SparseMatrix<double> &prm1,Eigen::SparseMatrix<double> &prm2){
		Eigen::SparseMatrix<double> temp( SumOfPoint ,  1 );
		for(int i=0;i< SumOfPoint ;i++)
			temp.coeffRef(i , 0) = prm1.coeffRef(point_num , i);
		temp.coeffRef(point_num , 0) = 0;
		temp = val  *  temp;
		prm2.coeffRef(point_num , 0) = val;
		prm2 = prm2 - temp;
		for(int i=0;i< SumOfPoint ;i++){
			if( prm1.coeffRef(i , point_num)==0 ) continue;
			prm1.coeffRef(i ,point_num)=0;
		}
		for(int i=0;i< SumOfPoint ;i++){
			if( prm1.coeffRef(point_num , i)==0 ) continue;
			prm1.coeffRef(point_num , i)=0;
		}
		prm1.coeffRef(point_num ,point_num) = 1;
		return;
	}

	void SetBC2D(
		int num
		,bool InverseFlag
		,Eigen::SparseMatrix<double> &prm1
		,Eigen::SparseMatrix<double> &prm2){
		//Generate Matrix
		prm2.resize(SumOfPoint,1);
		if(InverseFlag==false){	//ElectrodeSet Inverse Flag
			prm2.coeffRef(this->esp->Get(num,0),0) = 2.0;
			prm2.coeffRef(this->esp->Get(num,1),0) = -2.0;
		}else{
			prm2.coeffRef(this->esp->Get(num,2),0) = 2.0;
			prm2.coeffRef(this->esp->Get(num,3),0) = -2.0;
		}
		prm1 = FEM_base;
		SetDirchletCondition(
			this->esp->adp->ConvertPoint(
				this->esp->adp->NumOfEleX/2
				,this->esp->adp->NumOfEleY/2)
			,0.0,prm1,prm2);
		prm1.prune(0,0);
		return;
	}

	void SetBC3D(
		int num
		,bool InverseFlag
		,Eigen::SparseMatrix<double> &prm1
		,Eigen::SparseMatrix<double> &prm2){
		//Generate Matrix
		prm2.resize(SumOfPoint,1);
		if(InverseFlag==false){	//ElectrodeSet Inverse Flag
			prm2.coeffRef(this->esp->Get(num,0),0) = 2.0;
			prm2.coeffRef(this->esp->Get(num,1),0) = -2.0;
		}else{
			prm2.coeffRef(this->esp->Get(num,2),0) = 2.0;
			prm2.coeffRef(this->esp->Get(num,3),0) = -2.0;
		}
		prm1 = FEM_base;
		SetDirchletCondition(
			this->esp->adp->ConvertPoint(
				this->esp->adp->NumOfEleX/2
				,this->esp->adp->NumOfEleY/2
				,this->esp->adp->NumOfEleZ/2),
			0.0,prm1,prm2);
		prm1.prune(0,0);
		return;
	}

	Eigen::VectorXd Solve(Eigen::SparseMatrix<double> &prm1,Eigen::SparseMatrix<double> &prm2){
		Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > cg(prm1);
		Eigen::MatrixXd temp(prm2);
		Eigen::VectorXd ans = cg.solve(temp);
		/*
		std::cout << "#iterations:     " << cg.iterations() << std::endl;
		std::cout << "estimated error: " << cg.error()      << std::endl;
		*/
		return ans;
	}

	Eigen::VectorXd ImpedanceSense2D(Eigen::VectorXd &PotenA,Eigen::VectorXd &PotenB){
		Eigen::VectorXd res;
		for(int i=0;i<SumOfElement;i++){
			for(int j=0;j<4;j++){
					_ele[i].fai(j)= PotenA(_ele[i].point_num[j]);
					_ele[i].psi(j)= PotenB(_ele[i].point_num[j]);
			}
		}
		res = Eigen::VectorXd::Zero( SumOfElement );
		for(int id=0;id<SumOfElement;id++){
			for(int i=0;i<gi.bunten_num;i++){
				for(int j=0;j<gi.bunten_num;j++){
						double n1 = gi.bunten[i];
						double n2 = gi.bunten[j];
						double t_weight = gi.weight[i] * gi.weight[j];
						double t_res = 0.0;
						t_res += _ele[id].Dim2_CalcFaiDiff_fai(0,n1,n2) * _ele[id].Dim2_CalcFaiDiff_psi(0,n1,n2);
						t_res += _ele[id].Dim2_CalcFaiDiff_fai(1,n1,n2) * _ele[id].Dim2_CalcFaiDiff_psi(1,n1,n2);
						res(id) += t_res * t_weight;
				}
			}
		}
		return res;
	}

	Eigen::VectorXd ImpedanceSense3D(Eigen::VectorXd &PotenA,Eigen::VectorXd &PotenB){
		Eigen::VectorXd res;
		for(int i=0;i<SumOfElement;i++){
			for(int j=0;j<8;j++){
				_ele[i].fai(j)= PotenA(_ele[i].point_num[j]);
				_ele[i].psi(j)= PotenB(_ele[i].point_num[j]);
			}
		}
		res = Eigen::VectorXd::Zero( SumOfElement );
		for(int id=0;id<SumOfElement;id++){
			for(int i=0;i<gi.bunten_num;i++){
				for(int j=0;j<gi.bunten_num;j++){
					for(int k=0;k<gi.bunten_num;k++){
						double n1 = gi.bunten[i];
						double n2 = gi.bunten[j];
						double n3 = gi.bunten[k];
						double t_weight = gi.weight[i] * gi.weight[j] * gi.weight[k];
						double t_res = 0.0;
						t_res += _ele[id].Dim3_CalcFaiDiff_fai(0,n1,n2,n3) * _ele[id].Dim3_CalcFaiDiff_psi(0,n1,n2,n3);
						t_res += _ele[id].Dim3_CalcFaiDiff_fai(1,n1,n2,n3) * _ele[id].Dim3_CalcFaiDiff_psi(1,n1,n2,n3);
						t_res += _ele[id].Dim3_CalcFaiDiff_fai(2,n1,n2,n3) * _ele[id].Dim3_CalcFaiDiff_psi(2,n1,n2,n3);
						res(id) += t_res * t_weight;
					}
				}
			}
		}
		return res;
	}

	Eigen::VectorXd CalcZ(){	// For Translate-Filter
		std::cout << "**Making Kanzan Filter" << std::endl;
		Eigen::VectorXd SigmaMap = Eigen::VectorXd::Zero( SumOfElement );
		for(int i=0;i<SumOfElement;i++) SigmaMap(i) = 1.0;
		return CalcZ( SigmaMap );
	}

	Eigen::VectorXd CalcZ(Eigen::VectorXd &SigmaMap){
		Eigen::VectorXd res(esp->NumOfES);
		Eigen::SparseMatrix<double> FEM_lhs;
		Eigen::SparseMatrix<double> FEM_rhs;
		std::cout << "**Ready Assumed Sigma" << std::endl;
		if( this->esp->adp->NumOfDim == 2 ) SetSigma2D(SigmaMap);
		if( this->esp->adp->NumOfDim == 3 ) SetSigma3D(SigmaMap);
		std::cout << "**Calc" << std::endl;
		for(int i=0;i<(esp->NumOfES);i++){
			if( this->esp->adp->NumOfDim == 2 ) SetBC2D(i,false,FEM_lhs,FEM_rhs);
			if( this->esp->adp->NumOfDim == 3 ) SetBC3D(i,false,FEM_lhs,FEM_rhs);
			Eigen::VectorXd tmp = Solve(FEM_lhs,FEM_rhs);
			res(i) = tmp(esp->Get(i,3)) - tmp(esp->Get(i,2));
		}
		return res;
	}
	
	void CalcISM(Eigen::VectorXd &SigmaMap,Eigen::MatrixXd &ISM,Eigen::VectorXd &DeltaZ){
		Eigen::SparseMatrix<double> FEM_lhs;
		Eigen::SparseMatrix<double> FEM_rhs;
		Eigen::VectorXd PotenA;
		Eigen::VectorXd PotenB;
		std::cout << "**Ready Assumed Sigma" << std::endl;
		if( this->esp->adp->NumOfDim == 2 ) SetSigma2D(SigmaMap);
		if( this->esp->adp->NumOfDim == 3 ) SetSigma3D(SigmaMap);
		std::cout << "**Calc" << std::endl;
		for(int i=0;i<(esp->NumOfES);i++){
			if( this->esp->adp->NumOfDim == 2 ) SetBC2D(i,false,FEM_lhs,FEM_rhs);	//fai
			if( this->esp->adp->NumOfDim == 3 ) SetBC3D(i,false,FEM_lhs,FEM_rhs);	//fai
			PotenA = Solve(FEM_lhs,FEM_rhs);
			if( this->esp->adp->NumOfDim == 2 ) SetBC2D(i,true,FEM_lhs,FEM_rhs);	//psi
			if( this->esp->adp->NumOfDim == 3 ) SetBC3D(i,true,FEM_lhs,FEM_rhs);	//psi
			PotenB = Solve(FEM_lhs,FEM_rhs);
			//calc
			DeltaZ(i) -= PotenA(esp->Get(i,3));
			DeltaZ(i) += PotenA(esp->Get(i,2));
			if( this->esp->adp->NumOfDim == 2 )
			ISM.row(i) = ImpedanceSense2D(PotenA,PotenB);
			if( this->esp->adp->NumOfDim == 3 )
			ISM.row(i) = ImpedanceSense3D(PotenA,PotenB);
		}
		return;
	}

};
