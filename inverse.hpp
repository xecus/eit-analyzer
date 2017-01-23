
class InverseProblem{
private:
public:

	Fem *femp;
	Eigen::VectorXd SigmaMap_Assume;
	Eigen::VectorXd PhantomZ;
	Eigen::VectorXd DeltaZ;
	Eigen::MatrixXd ISM;
	Eigen::MatrixXd DeltaSigma;

	~InverseProblem(){
	}
	InverseProblem(Fem *_femp){
		this->femp = _femp;
	}

	void StableSolve(){
		//std::cout << "**SV Decomposing..." << std::endl;
		int ism_rows = this->ISM.rows();
		int ism_cols = this->ISM.cols();
		//std::cout << "ISM=(" << ism_rows << "," << ism_cols << ")" << std::endl;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(this->ISM, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::MatrixXd S = svd.singularValues();
		Eigen::MatrixXd S2 = Eigen::MatrixXd::Zero(ism_rows,ism_cols);
		for(int i=0;i<ism_cols;i++) S2(i,i) = S(i,0);
		Eigen::MatrixXd U = svd.matrixU();
		Eigen::MatrixXd V = svd.matrixV();
		/*
		double lambda_temp = 0.0;
		Eigen::FullPivLU< Eigen::MatrixXd > lu( S2.transpose() * S2
			+ lambda_temp * lambda_temp * Eigen::MatrixXd::Identity(ism_cols,ism_cols)  );
		Eigen::MatrixXd PseudoInverse = V * lu.inverse() * S2.transpose() * U.transpose();
		Eigen::VectorXd ans = PseudoInverse * this->DeltaZ;
		Eigen::VectorXd err = this->ISM * ans - this->DeltaZ;
		this->DeltaSigma = ans;
		*/
		//std::cout << "**Try to changing lambda..." << std::endl;
		double lambda_start = -5;
		double lambda_stop = 5;
		double lambda_step = sqrt(10);
		int TryCnt=0;
		for(double lambda_temp = pow(10.0,(double)lambda_start);
			lambda_temp < pow(10.0,(double)lambda_stop);
			lambda_temp *= lambda_step
			){
			//std::cout << "**lambda[" << lambda_temp << "]..." << std::endl;
			Eigen::FullPivLU< Eigen::MatrixXd > lu( S2.transpose() * S2
				+ lambda_temp * lambda_temp * Eigen::MatrixXd::Identity(ism_cols,ism_cols)  );
			Eigen::MatrixXd PseudoInverse = V * lu.inverse() * S2.transpose() * U.transpose();
			Eigen::VectorXd ans = PseudoInverse * this->DeltaZ;
			Eigen::VectorXd err = this->ISM * ans - this->DeltaZ;
			this->DeltaSigma = ans;

			double th = femp->SumOfElement * (0.5 * 0.5);
			if(ans.norm() < th){
				if(TryCnt==0 && lambda_start > -10)
					lambda_start-=2;
				/*
				std::cout << "info:Regulation Parameter=" << lambda_temp << std::endl;
				std::cout << "info:ans_temp.norm()=" << ans_temp.norm() << std::endl;
				std::cout << "info:err_temp.norm()=" << err_temp.norm() << std::endl;
				*/
				break;
			}

			TryCnt++;
		}

		return;
	}

	InverseProblem *DoInverse(){
		//Make ISM
		//std::cout << "**Copy PhantomZ -> DeltaZ" << std::endl;
		DeltaZ = PhantomZ;
		//std::cout << "**Making ISM" << std::endl;
		femp->CalcISM(SigmaMap_Assume,ISM,DeltaZ);
		//std::cout << "**Call Stable Sovler" << std::endl;
		StableSolve();
		//Update
		//std::cout << "**Updating Sigma Distribution" << std::endl;
		SigmaMap_Assume += DeltaSigma;
		for(int i=0;i<(this->femp->SumOfElement);i++){
			if(SigmaMap_Assume(i) >= 5.0) SigmaMap_Assume(i) = 5.0;
			if(SigmaMap_Assume(i) <= 0.5) SigmaMap_Assume(i) = 0.5;
		}
		//std::cout << "[SigmaMap] ";
		//for(int i=0;i<5;i++) std::cout << SigmaMap_Assume(i) << " ";
		//std::cout << std::endl;
		return this;
	}

	InverseProblem *LoadPhantomZ(){
		
		Eigen::VectorXd SigmaMap_True;

		//std::cout << "**Define True Map" << std::endl;
		SigmaMap_True = Eigen::VectorXd::Zero( this->femp->SumOfElement );
		for(int i=0;i<(this->femp->SumOfElement);i++) SigmaMap_True(i) = 1.5;
		SigmaMap_True(femp->esp->adp->ConvertElement(1,1)) = 2.0;
		SigmaMap_True(femp->esp->adp->ConvertElement(2,1)) = 2.0;
		SigmaMap_True(femp->esp->adp->ConvertElement(1,2)) = 2.0;
		SigmaMap_True(femp->esp->adp->ConvertElement(2,2)) = 2.0;
		SigmaMap_True(femp->esp->adp->ConvertElement(5,5)) = 3.0;
		//std::cout << "**Calc PhantomZ" << std::endl;
		PhantomZ = femp->CalcZ(SigmaMap_True);
		return this;
	}

	InverseProblem *Init(){
		//std::cout << "**Allocating Eigen Matrix" << std::endl;
		SigmaMap_Assume = Eigen::VectorXd::Zero( this->femp->SumOfElement );
		PhantomZ = Eigen::VectorXd::Zero( femp->esp->NumOfES );
		DeltaZ = Eigen::VectorXd::Zero( femp->esp->NumOfES );
		ISM = Eigen::MatrixXd::Zero( femp->esp->NumOfES , femp->SumOfElement );
		DeltaSigma = Eigen::VectorXd::Zero( this->femp->SumOfElement );
		for(int i=0;i<(this->femp->SumOfElement);i++) SigmaMap_Assume(i) = 1.0;
		//std::cout << "Unknown=" << femp->SumOfElement << std::endl;
		//std::cout << "NumOfES=" << femp->esp->NumOfES << std::endl;
		return this;
	}

	/*
	void SavePotential(std::string fn,Eigen::VectorXd &prm1){

		std::ofstream ofs;
		ofs.open( fn.c_str() );

		ofs << "[Debug Mode]" << std::endl;
		ofs << "NumOfDim" << "," << this->esp->adp->NumOfDim << std::endl; 
		ofs << "NumOfEleX" << "," << this->esp->adp->NumOfEleX << std::endl;
		ofs << "NumOfEleY" << "," << this->esp->adp->NumOfEleY << std::endl;
		if(this->esp->adp->NumOfDim==3)
		ofs << "NumOfEleZ" << "," << this->esp->adp->NumOfEleZ << std::endl;
		ofs << "NumOfPntX" << "," << this->esp->adp->NumOfEleX+1 << std::endl; 
		ofs << "NumOfPntY" << "," << this->esp->adp->NumOfEleY+1 << std::endl; 
		if(this->esp->adp->NumOfDim==3)
		ofs << "NumOfPntZ" << "," << this->esp->adp->NumOfEleZ+1 << std::endl; 

		if(this->esp->adp->NumOfDim==2){
			int x = this->esp->adp->NumOfEleX+1;
			int y = this->esp->adp->NumOfEleY+1;
			ofs << std::endl;
			for(int j=0;j<y;j++){
				for(int k=0;k<x;k++){
					int id = this->esp->adp->ConvertPoint(k,j);
					ofs << prm1(id);
					if(k==x-1) ofs << std::endl; else ofs << ",";
				}
			}
			ofs.close();
		}
		if(this->esp->adp->NumOfDim==3){
			int x = this->esp->adp->NumOfEleX+1;
			int y = this->esp->adp->NumOfEleY+1;
			int z = this->esp->adp->NumOfEleZ+1;
			ofs << std::endl;
			for(int i=0;i<z;i++){
				ofs << "[" << i << "]" << std::endl;
				for(int j=0;j<y;j++){
					for(int k=0;k<x;k++){
						int id = this->esp->adp->ConvertPoint(k,j,i);
						ofs << prm1(id);
						if(k==x-1) ofs << std::endl; else ofs << ",";
					}
				}
				ofs << std::endl;
			}
			ofs.close();
		}

		return;
	}
	*/

};
