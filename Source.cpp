#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include<random>

using namespace std;
using namespace matlab::engine;
using boost::math::normal;

double sum(vector<double>s) {
	/*This function returns sum of vector elements */
	double out = 0;
	for (int i=0; i < s.size(); i++) {
		out = out + s[i];
	}
	return out;
}

vector<double> dot(vector<double>x, vector<double>y) {
	/*This function returns the elemntwise product of two vector */
	vector<double>out;
	for (int i=0; i < x.size(); i++) {
		out.push_back(x[i] * y[i]);
	}
	return out;
}


vector<vector<double>> genrate() {
	vector<vector<double>> out;    /*Container to store the state and the measurement signals*/
	vector<double>x;               /*state signal*/
	vector<double>y;               /*measurement signal*/
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*The following lines of code are for implementing the noise generator for the state and measurement equation*/
	default_random_engine generator{};
	normal_distribution<double> distribution1(0, sqrt(10));
	normal_distribution<double> distribution2(0, 1);
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*Evaluting the state and the measurement signals*/
	double x0 = distribution2(generator);
	for (int i = 0; i < 100; i++) {
		double temp = .5*x0 + (25 * x0) / (1 + pow(x0, 2)) + distribution1(generator) + 8 * cos(1.2*i);
		x.push_back(temp);
		x0 = x[i];
		double temp1 = pow(x[i], 2) / 20 + distribution2(generator);
		y.push_back(temp1);
	}
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*Putting "x" and "y " in the container "out"*/
	out.push_back(x);
	out.push_back(y);
	return out;
}








vector<vector<double>> ourfun(vector<double>y) {
	normal pd1(0, sqrt(10)); /*Pdf of the state equation noise*/
	normal pd2;              /*Pdf of measurement equation noise*/
	vector<double>x01;       /*Grid values vector*/
	vector<double> w0;       /*Temporary container*/
	vector<double> wprior;   /*Prior estimated pdf*/
	vector<double> g;        /*Temporary container*/
	vector<double> wpost;    /*Post estimated pdf*/
	vector<double> tempy;    /*Temporary container*/
	vector<vector<double>>out; /*Container to store the output */


	/* Create the grid values vector "x01" */
	for (double i = -25; i <= 25; i = i + .5) {
		x01.push_back(i);
	}
	
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*Evaluate the probability density function at each grid value , and store it in the "w0"*/
	for (int h = 0; h < x01.size(); h++) {
		w0.push_back(pdf(pd1, x01[h]));
	}
	
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*Outer for loop to estimate 100 points of the state*/
	for (int outer = 0; outer < 100; outer++) {
		
		/*for loop to evalute the prior estimated pdf*/
		for (int i1 = 0; i1 < x01.size(); i1++) {
			g.clear();

			for (int i2 = 0; i2 < x01.size(); i2++) {
				double temp = x01[i1] - (x01[i2] / 2 + (25 * x01[i2]) / (1 + pow(x01[i2], 2)) + 8 * cos(1.2*(outer+1)));
				g.push_back(temp);
			}
			for (int i3 = 0; i3 < g.size(); i3++) {
				g[i3] = pdf(pd1, g[i3]);

			}
			vector<double> temp = dot(w0, g);
			wprior.push_back(sum(temp));

		}
/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
      /*The following lines of code are for evaluating the posterior pdf */


		tempy.clear();
		for (int o = 0; o < x01.size(); o++) {
			double temp = y[outer] - pow(x01[o], 2) / 20;
			tempy.push_back(temp);

	}
		for (int o1 = 0; o1 < tempy.size(); o1++) {
			tempy[o1] = pdf(pd2, tempy[o1]);
		}
	
		wpost = dot(tempy, wprior);
		double ss = sum(wpost);
		for (int k = 0; k < wpost.size(); k++) {
			wpost[k] = wpost[k] / ss;
		}
		
		/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		
		/*Put each estimated Posterior pdf in the "out" container*/
	
		out.push_back(wpost);
	
		/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
		/*Do the required preparations for the next iteration of the Outer for loop*/
	
		w0 = wpost;
	
		wprior.clear();
	
	}/*End of the Outer for loop*/

	out.push_back(x01); /*Put "x01" in the final row of the cointainer "out"*/

	return out;






}














void main() {
	vector<vector<double>>input;    /*Container to store the signals that will  be passed to the filter*/
	vector<vector<double>>pdf_est; /*Container to store the output from the filter*/
	vector<double>xT;              /*True values for the estimated state*/
	vector<double>y;               /*Measurement signal*/
	vector<double>x01;             /*Grid values vector*/
	vector<double>esti;            /* The estimated state */
	
	
	
	input = genrate();             /*genrate function returns the true state(xT) and the measurement signal(y)*/
	xT = input[0];
	y = input[1];
	pdf_est = ourfun(y);          /*ourfun function is the grid filter*/
	x01 = pdf_est[pdf_est.size() - 1];
	
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
	/*In here , we evalute the mean value for each resulted pdf and store the value in "esti" vector */
	for (int i = 0; i < pdf_est.size() - 1; i++) {
		vector<double>temp = pdf_est[i];
		double temp1 = sum(dot(x01, temp));
		esti.push_back(temp1);
	}
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
	/*The next lines of code are responsible  for the following:
1-	Create a MATALB engine instance .
2-	Put each of “xT” and “esti ” in arrays supported by MATLAB. These arrays are “xt_mat” & “esti_mat” respectively.
3-	Send these arrays to MATLAB engine to plot them. 
 */
	unique_ptr<MATLABEngine> matlabPtr=startMATLAB(); 
	matlab::data::ArrayFactory factory;               
	auto est_mat = factory.createArray<double>({ 1, esti.size() });
	auto xt_mat= factory.createArray<double>({ 1, xT.size() });
	for (int i = 0; i < xT.size(); i++) {
		est_mat[i] = esti[i];
		xt_mat[i] = xT[i];

	}
	matlabPtr->setVariable(u"est", move(est_mat));
	matlabPtr->setVariable(u"xtt", move(xt_mat));
	matlabPtr->eval(u"plot([1:100],xtt,'b',[1:100],est,'r')");

	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
	system("pause");

}