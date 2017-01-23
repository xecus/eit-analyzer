
class AnalysisDomain{
private:
public:
	//Domain
	int NumOfDim;
	int NumOfEleX;
	int NumOfEleY;
	int NumOfEleZ;
	double LenghOfX;
	double LenghOfY;
	double LenghOfZ;
	//ES
	int ModeOfES;
	int PosOfE_ZDir;
	int NumOfE_ZDir;
	~AnalysisDomain(){
	}
	AnalysisDomain(){
	}

	//Experiment Model
	AnalysisDomain *Test2D_Exp1(){
		//std::cout << "**Loading Phantom Experiment Model(20150326)" << std::endl;
		this->NumOfDim = 2;
		this->NumOfEleX = 10;
		this->NumOfEleY = 10;
		this->NumOfEleZ = 0;	//dont care
		this->LenghOfX = 0.15;
		this->LenghOfY = 0.15;
		this->LenghOfZ = 0.0;	//dont care
		this->ModeOfES = 1;
		this->PosOfE_ZDir = 0;	//fix
		this->NumOfE_ZDir = 1;	//fix
		return this;
	}

	//Test For 2D
	AnalysisDomain *Test2D(){
		this->NumOfDim = 2;
		this->NumOfEleX = 8;
		this->NumOfEleY = 8;
		this->NumOfEleZ = 0;	//dont care
		this->LenghOfX = 0.15;
		this->LenghOfY = 0.15;
		this->LenghOfZ = 0.0;	//dont care
		this->ModeOfES = 0;
		this->PosOfE_ZDir = 0;	//fix
		this->NumOfE_ZDir = 1;	//fix
		return this;
	}

	//Test For 3D
	AnalysisDomain *Test3D(){
		this->NumOfDim = 3;
		this->NumOfEleX = 8;
		this->NumOfEleY = 8;
		this->NumOfEleZ = 8;
		this->LenghOfX = 0.15;
		this->LenghOfY = 0.15;
		this->LenghOfZ = 0.08;
		this->ModeOfES = 0;
		this->PosOfE_ZDir = this->NumOfEleZ/2 - 2;
		this->NumOfE_ZDir = 3;
		return this;
	}

	//User
	AnalysisDomain *Set(int _NumOfDim,int _NumOfEle[],double _LengthOf[],int _ES[]){
		this->NumOfDim = _NumOfDim;
		this->NumOfEleX = _NumOfEle[0];
		this->NumOfEleY = _NumOfEle[1];
		this->NumOfEleZ = _NumOfEle[2];
		this->LenghOfX = _LengthOf[0];
		this->LenghOfY = _LengthOf[1];
		this->LenghOfZ = _LengthOf[2];
		this->ModeOfES = _ES[0];
		this->PosOfE_ZDir = _ES[1];
		this->NumOfE_ZDir = _ES[2];
		if( _NumOfDim==2 ){
			this->PosOfE_ZDir = 0;
			this->NumOfE_ZDir = 1;
		}
		return this;
	}


	//Converter(Utilitiy)
	int ConvertPoint(int x,int y,int z){
		int res = ( z * ((this->NumOfEleX+1)*(this->NumOfEleY+1)) + y *(this->NumOfEleX+1) + x );
		if( x > this->NumOfEleX || y > this->NumOfEleY || z > this->NumOfEleZ ){
			//std::cout << "(ConvertPoint) Err" << std::endl;
			res = -1;
		}
		return res;
	}
    int ConvertElement(int x,int y,int z){
		int res = ( z * (this->NumOfEleX*this->NumOfEleY) + y * this->NumOfEleX + x );
		if( x >= this->NumOfEleX || y >= this->NumOfEleY || z >= this->NumOfEleZ ){
			//std::cout << "(ConvertElement) Err" << std::endl;
			res = -1;
		}
		return res;
    }
    int ConvertPoint(int x,int y){
		int res = ( y * (this->NumOfEleX+1) + x );
		if( x > this->NumOfEleX || y > this->NumOfEleY ){
			//std::cout << "(ConvertPoint) Err" << std::endl;
			res = -1;
		}
		return res;
    }
    int ConvertElement(int x,int y){
		int res = ( y *NumOfEleX + x );
		if( x >= this->NumOfEleX || y >= this->NumOfEleY ){
			//std::cout << "(ConvertElement) Err" << std::endl;
			res = -1;
		}
		return res;
    }

};
