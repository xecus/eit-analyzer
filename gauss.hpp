
class GaussInt{
private:
public:
	int bunten_num;
	double *bunten;
	double *weight;
	void Gauss3(){
		bunten_num = 3;
		bunten = new double[bunten_num];
		weight = new double[bunten_num];
		bunten[0]=-0.7745966692414833770359;
		weight[0]=+0.555555555555555555556;
		bunten[1]=+0.0;
		weight[1]=+0.8888888888888888888889;
		bunten[2]=+0.7745966692414833770359;
		weight[2]=+0.555555555555555555556;
		return;
	}
	void Gauss10(){
		bunten_num = 10;
		bunten = new double[bunten_num];
		weight = new double[bunten_num];
		bunten[0]=-0.973906528517171720078;
		weight[0]=0.0666713443086881375936;
		bunten[1]=-0.8650633666889845107321;
		weight[1]=0.1494513491505805931458;
		bunten[2]=-0.6794095682990244062343;
		weight[2]=0.219086362515982043996;
		bunten[3]=-0.4333953941292471907993;
		weight[3]=0.269266719309996355091;
		bunten[4]=-0.1488743389816312108848;
		weight[4]=0.295524224714752870174;
		bunten[5]=0.1488743389816312108848;
		weight[5]=0.295524224714752870174;
		bunten[6]=0.4333953941292471907993;
		weight[6]=0.269266719309996355091;
		bunten[7]=0.6794095682990244062343;
		weight[7]=0.219086362515982043996;
		bunten[8]=0.8650633666889845107321;
		weight[8]=0.1494513491505805931458;
		bunten[9]=0.973906528517171720078;
		weight[9]=0.0666713443086881375936;
		return;
	}
	~GaussInt(){
		delete []bunten;
		delete []weight;
		return;
	}
	GaussInt(){
		Gauss10();
	}
};
GaussInt gi;
