
class ESet{

private:
public:

	AnalysisDomain *adp;

	int NumOfE_Layer;
	int NumOfE_Sum;
	int NumOfES;
	int *EsArray;

	~ESet(){
	}
	ESet(){
	}

	//Print Detail
	void Print(){
		//std::cout << "[Eset Print]" << std::endl;
		//std::cout << "NumOfE_Layer=" << NumOfE_Layer << std::endl;
		//std::cout << "NumOfE_Sum=" << NumOfE_Sum << std::endl;
		//std::cout << "NumOfES=" << NumOfES << std::endl;
		return;
	}


	ESet(AnalysisDomain *_adp){
		this->adp = _adp;
		//std::cout << "**Calc Electrode Set Index" << std::endl;
		//Electrode Set Type I
		if(adp->ModeOfES == 0)
			this->NumOfE_Layer = (adp->NumOfEleX/2 + adp->NumOfEleY/2)*2;
		//Electrode Set Type II
		if(adp->ModeOfES == 1)
			this->NumOfE_Layer = (adp->NumOfEleX-1 + adp->NumOfEleY-1)*2;
		//Electrode Set Type III
		if(adp->ModeOfES == 2)
			this->NumOfE_Layer = (adp->NumOfEleX/2-1 + adp->NumOfEleY/2-1)*2;
		this->NumOfE_Sum = this->NumOfE_Layer * adp->NumOfE_ZDir;
		this->NumOfES = (this->NumOfE_Sum-1)*(this->NumOfE_Sum-2)/2;
		EsArray = new int[ this->NumOfES * 4 ];
		/*
		std::cout << "this->NumOfE_Layer=" << this->NumOfE_Layer << std::endl;
		std::cout << "this->NumOfE_Sum=" << this->NumOfE_Sum << std::endl;
		std::cout << "this->NumOfES=" << this->NumOfES << std::endl;
		*/
		//std::cout << "**Generating Electrode Set Array" << std::endl;
		Generate();
		//std::cout << "**Saving Electrode Set" << std::endl;
		SaveES("./out/ES.csv");
		return;
	}
	ESet(AnalysisDomain *_adp,std::string fn,int _numofes){
		this->adp = _adp;
		//std::cout << "**Calc Electrode Set Index" << std::endl;
		//Electrode Set Type I
		if(adp->ModeOfES == 0)
			this->NumOfE_Layer = (adp->NumOfEleX/2 + adp->NumOfEleY/2)*2;
		//Electrode Set Type II
		if(adp->ModeOfES == 1)
			this->NumOfE_Layer = (adp->NumOfEleX-1 + adp->NumOfEleY-1)*2;
		//Electrode Set Type III
		if(adp->ModeOfES == 2)
			this->NumOfE_Layer = (adp->NumOfEleX/2-1 + adp->NumOfEleY/2-1)*2;
		this->NumOfE_Sum = this->NumOfE_Layer * adp->NumOfE_ZDir;
		this->NumOfES = _numofes;
		EsArray = new int[ _numofes * 4 ];
		/*
		std::cout << "this->NumOfE_Layer=" << this->NumOfE_Layer << std::endl;
		std::cout << "this->NumOfE_Sum=" << this->NumOfE_Sum << std::endl;
		std::cout << "this->NumOfES=" << this->NumOfES << std::endl;
		*/

		int cnt = 0;
		int cnt2 = 0;
		int flag = 0;
		std::ifstream ifs(fn.c_str());
		std::string str;
		while( getline( ifs, str ) ){
			std::string token;
			std::stringstream stream( str );
			while( getline( stream, token, '\t' ) ) {
				if(token=="") continue;
				if(cnt%4==0){
					if(flag) flag=0; else flag=1;
					//cnt2=0;
				}
				if(flag==0){
					std::stringstream ss;
					ss << token;
					int tmp;
					ss >> tmp;
					//std::cout << "[" << cnt2 << "]token=" << tmp << std::endl;
					EsArray[cnt2++]=tmp;
				}
				cnt++;
			}

		}

		for(int i=0;i<_numofes;i++){
			//std::cout << "[" << i << "] ";
			//std::cout << EsArray[4*i+0] << " ";
			//std::cout << EsArray[4*i+1] << " ";
			//std::cout << EsArray[4*i+2] << " ";
			//std::cout << EsArray[4*i+3] << " ";
			//std::cout << std::endl;
		}

		//std::cout << "**Generating Electrode Set Array" << std::endl;
		Generate();
		//std::cout << "**Saving Electrode Set" << std::endl;
		SaveES("./out/ES.csv");
		return;
	}

	/* Electrode Number to Point Number */
	int ESNum2PntNum(int num){
		if(adp->ModeOfES == 0)
			return ESNum2PntNum_T1(num);
		if(adp->ModeOfES == 1)
			return ESNum2PntNum_T2(num);
		if(adp->ModeOfES == 2)
			return ESNum2PntNum_T3(num);
		return -1;
	}
	int ESNum2PntNum_T1(int num){
		int temp[4];
		temp[0] = adp->NumOfEleX/2;
		temp[1] = temp[0] + adp->NumOfEleY/2;
		temp[2] = temp[1] + adp->NumOfEleX/2;
		temp[3] = temp[2] + adp->NumOfEleY/2;
		int z = num / this->NumOfE_Layer + adp->PosOfE_ZDir;
		num = num % this->NumOfE_Layer;
		if(num < temp[0] ){
				return adp->ConvertPoint(1+num*2,0,z);
		}else if(num < temp[1] ){
			num -= temp[0];
			return adp->ConvertPoint(adp->NumOfEleX,1+num*2,z);
		}else if(num < temp[2] ){
			num -= temp[1];
			return adp->ConvertPoint(adp->NumOfEleX-(1+num*2),adp->NumOfEleY,z);
		}else if(num < temp[3] ){
			num -= temp[2];
			return adp->ConvertPoint(0,adp->NumOfEleY-(1+num*2),z);
		}else{
			return -1;
		}
	}
	int ESNum2PntNum_T2(int num){
		int temp[4];
		temp[0] = adp->NumOfEleX - 1;
		temp[1] = temp[0] + (adp->NumOfEleY-1);
		temp[2] = temp[1] + (adp->NumOfEleX-1);
		temp[3] = temp[2] + (adp->NumOfEleY-1);
		int z = num / this->NumOfE_Layer + adp->PosOfE_ZDir;
		num = num % this->NumOfE_Layer;
		if(num < temp[0] ){
			return adp->ConvertPoint(1+num,0,z);
		}else if(num < temp[1] ){
			num -= temp[0];
			return adp->ConvertPoint(adp->NumOfEleX,1+num,z);
		}else if(num < temp[2] ){
			num -= temp[1];
			return adp->ConvertPoint(adp->NumOfEleX-(1+num),adp->NumOfEleY,z);
		}else if(num < temp[3] ){
			num -= temp[2];
			return adp->ConvertPoint(0,adp->NumOfEleY-(1+num),z);
		}else{
			return -1;
		}
	}
	int ESNum2PntNum_T3(int num){
		int temp[4];
		temp[0] = (adp->NumOfEleX/2-1);
		temp[1] = temp[0] + (adp->NumOfEleY/2-1);
		temp[2] = temp[1] + (adp->NumOfEleX/2-1);
		temp[3] = temp[2] + (adp->NumOfEleY/2-1);
		int z = num / this->NumOfE_Layer + adp->PosOfE_ZDir;
		num = num % this->NumOfE_Layer;
		if(num < temp[0] ){
				return adp->ConvertPoint(2+num*2,0,z);
		}else if(num < temp[1] ){
			num -= temp[0];
			return adp->ConvertPoint(adp->NumOfEleX,2+num*2,z);
		}else if(num < temp[2] ){
			num -= temp[1];
			return adp->ConvertPoint(adp->NumOfEleX-(2+num*2),adp->NumOfEleY,z);
		}else if(num < temp[3] ){
			num -= temp[2];
			return adp->ConvertPoint(0,adp->NumOfEleY-(2+num*2),z);
		}else{
			return -1;
		}
	}

	/* Generate */
	int Generate(){
		int *tmp_set;
		tmp_set = new int[(this->NumOfE_Sum-1)*2];
		for(int i=0;i<this->NumOfE_Sum-1;i++){
			tmp_set[i*2] = i;
			tmp_set[i*2+1] = i+1;
		}
		int tempcnt=0;
		for(int i=0;i<this->NumOfE_Sum-1;i++){
			for(int j=i+1;j<this->NumOfE_Sum-1;j++){
				EsArray[tempcnt*4+0] = ESNum2PntNum(tmp_set[i*2]);
				EsArray[tempcnt*4+1] = ESNum2PntNum(tmp_set[i*2+1]);
				EsArray[tempcnt*4+2] = ESNum2PntNum(tmp_set[j*2]);
				EsArray[tempcnt*4+3] = ESNum2PntNum(tmp_set[j*2+1]);
				tempcnt++;
			}
		}
		delete[] tmp_set;
		return tempcnt;
	}

	int Get(int esNum,int info){
		if( esNum >= NumOfES ){
			//std::cout << "(Error!)" << std::endl;
			return -1;
		}else{
			return EsArray[ esNum *4 + info];
		}
	}

	ESet *SaveES(std::string fn){
		int m = NumOfES;
		int n = 4;
		std::ofstream ofs;
		ofs.open( fn.c_str() );
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				ofs << Get(i,j);
				if(j==n-1) ofs << std::endl;
				else ofs << ",";
			}
		}
		ofs.close();
		return this;
	}


};

