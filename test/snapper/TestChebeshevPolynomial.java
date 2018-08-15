package snapper;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.TestCase;

public class TestChebeshevPolynomial extends TestCase {
	
	@Test
	public void testChebeshevPolynomial() {
		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(60, 100);
		
		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		
		// convert factors to function values
		p.aToF();
		
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < f.length; i++) {
			if (Math.abs(f[i] - p.f[i]) > maxErr) {
				maxErr = Math.abs(f[i] - p.f[i]);
				iErr = i;
			}
		}
		System.out.println("Max error: " + maxErr + "\nat entry " + iErr);
		assertTrue(Math.abs(maxErr) < 1e-14);
	}

	
	@Test
	public void testChebeshevPolynomial2() {

		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(10, 32);

		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		System.out.println(Arrays.toString(f));
		double [] fin = new double[] {				
				0.00000000e+00,4.00424050e-19,3.49498706e-13,8.55254801e-10
				,1.75230607e-07,8.70224070e-06,1.67902916e-04,1.61580650e-03
				,8.98484427e-03,3.16714850e-02,7.52116961e-02,1.25324508e-01
				,1.50543597e-01,1.32619063e-01,8.64907774e-02,4.18856963e-02
				,1.50204264e-02,3.95490195e-03,7.53575618e-04,1.01716090e-04
				,9.44244301e-06,5.79210760e-07,2.22457764e-08,4.97232577e-10
				,5.84647512e-12,3.13273808e-14,6.18871196e-17,3.23947049e-20
				,2.57924641e-24,1.10157019e-29,2.42664019e-37,1.56377455e-50
				,0.00000000e+00				
				};
		for (int i = 0; i < N; i++) {
			System.out.println(fin[i] - f[i]);
		}
		double [] a = new double [] {
				2.10739037e-02,1.53264754e-02,-2.88498361e-02,-3.33576230e-02
				,1.80311475e-03,2.73147519e-02,1.59228114e-02,-9.11155969e-03
				,-1.58167632e-02,-3.70276402e-03,6.92228525e-03,5.84414927e-03
				,-2.73424845e-04,-2.87309624e-03,-1.38012637e-03,4.56719059e-04
				,7.45906917e-04,2.06888613e-04,-1.39383680e-04,-1.26441912e-04
				,-2.08314301e-05,2.15809534e-05,1.42207415e-05,1.88922381e-06
				,-1.72975078e-06,-1.02661462e-06,-1.86219733e-07,4.60513342e-08
				,3.70425070e-08,1.07994016e-08,1.79057434e-09,1.67866354e-10
				,6.99444148e-12	
		};
		System.out.println(Arrays.toString(p.a));
		for (int i = 0; i < N; i++) {
			System.out.println(a[i] - p.a[i]);
		}
		
		p.aToF();
		System.out.println(Arrays.toString(f));
		double [] fout = new double[]{
				
				-3.49721033e-12,3.49721422e-12,-3.14772324e-12,8.58752026e-10
				,1.75227109e-07,8.70224420e-06,1.67902913e-04,1.61580650e-03
				,8.98484426e-03,3.16714850e-02,7.52116961e-02,1.25324508e-01
				,1.50543597e-01,1.32619063e-01,8.64907774e-02,4.18856963e-02
				,1.50204264e-02,3.95490196e-03,7.53575614e-04,1.01716093e-04
				,9.44243951e-06,5.79214258e-07,2.22422792e-08,5.00729789e-10
				,2.34926488e-12,3.52854773e-12,-3.49715959e-12,3.49722334e-12
				,-3.49722106e-12,3.49723254e-12,-3.49722524e-12,3.49722640e-12
				,-3.49721727e-12
		};
		
		for (int i = 0; i < N; i++) {
			System.out.println(fout[i] - p.f[i]);
			assertEquals(fout[i], p.f[i], 5e-10);
		}
/* Python code for generating the above numbers in this test
 
 %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps

def spectral_to_nodes(c,N):
    c = np.reshape(c,(len(c),1)) #This is to get a column vector, since the fourier transform output depends on the orientation
    cc = np.append(c[1:],np.flip(c[1:-1], axis = 0))/2 #Need to make some changes in order to get the Chebyshev coefficients
    ap = np.append(c[0], cc) #Apply the FFT and throw away the complex values
    F=np.real(np.fft.fft(ap))
    function_values=F[0:N+1] 
    return function_values

def nodes_to_spectral(f,N):
    f1 = np.reshape(f,(len(f),1)) #Need input as a column vector
    f2 = np.reshape(np.flip(f[1:-1], axis = 0),(len(f)-2,1)) #Need to make some alteration before we apply the inverse FFT 
    f = np.real(np.fft.ifft(np.append(f1,f2))) #Apply the inverse fourier transform and throw away the complex values 
    spectral_coefficients = np.append(f[0], 2*f[1:N]) #Get the Chebyshev coefficients
    spectral_coefficients = np.append(spectral_coefficients, f[N])
    return spectral_coefficients

#Example output
n = 32
m = 10
N = 32
circle = -np.arange(0,np.pi+np.pi/n,np.pi/n)
print("circle=", circle)
x = 1-(np.cos(circle) + 1)/2
print("x=", x)
f = sps.binom.pmf(m, n, x, loc=0)
print("Node values in")
print(f)
c = nodes_to_spectral(f,N)
f = spectral_to_nodes(c,n)
print("Chebyshev coefficients")
print(c)
print("Node values out")
print(f)
plt.plot(x,f)
plt.show()
		
 */
	}
	
	@Test
	public void testChebeshevPolynomial3() {

		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(5, 20);

		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		System.out.println(Arrays.toString(f));
		double [] fin = new double[] {
				0.00000000e+00,1.20974891e-09,1.09793114e-06,5.17440974e-05
				,6.91861550e-04,4.46399871e-03,1.75885770e-02,4.79228431e-02
				,9.71182265e-02,1.53201567e-01,1.93733134e-01,2.00188845e-01
				,1.71090106e-01,1.21772116e-01,7.23753262e-02,3.58898848e-02
				,1.47857666e-02,5.02184976e-03,1.39012781e-03,3.08698324e-04
				,5.38364771e-05,7.16869967e-06,7.02003144e-07,4.80661306e-08
				,2.14677751e-09,5.66876154e-11,7.65624338e-13,4.22110842e-15
				,6.50511413e-18,1.37654800e-21,8.10105561e-27,8.11098433e-36
				,0.00000000e+00
				};
		for (int i = 0; i < N; i++) {
			System.out.println(fin[i] - f[i]);
		}
		double [] a = new double [] {
						3.55517978e-02,3.38588550e-02,-3.38588550e-02,-5.74128411e-02
						,-2.24498930e-02,2.10808176e-02,3.00313323e-02,1.05338660e-02
						,-7.49147614e-03,-9.88861313e-03,-3.59581830e-03,1.18897995e-03
						,1.86680569e-03,8.07130127e-04,3.55340308e-05,-1.41120865e-04
						,-8.43228190e-05,-2.65095150e-05,-5.07629011e-06,-5.64032234e-07
						,-2.82016117e-08,-6.36968776e-18,2.09251019e-17,1.17093835e-17
						,-7.80625564e-18,-1.40946282e-17,0.00000000e+00,1.56125113e-17
						,1.47451495e-17,-6.93889390e-18,-1.04083409e-17,-1.73472348e-18
						,-3.46944695e-18
		};
		System.out.println(Arrays.toString(p.a));
		for (int i = 0; i < N; i++) {
			System.out.println(a[i] - p.a[i]);
		}
		
		p.aToF();
		System.out.println(Arrays.toString(f));
		double [] fout = new double[]{
				-2.42861287e-17,1.20974890e-09,1.09793114e-06,5.17440974e-05
				,6.91861550e-04,4.46399871e-03,1.75885770e-02,4.79228431e-02
				,9.71182265e-02,1.53201567e-01,1.93733134e-01,2.00188845e-01
				,1.71090106e-01,1.21772116e-01,7.23753262e-02,3.58898848e-02
				,1.47857666e-02,5.02184976e-03,1.39012781e-03,3.08698324e-04
				,5.38364771e-05,7.16869967e-06,7.02003144e-07,4.80661306e-08
				,2.14677749e-09,5.66876060e-11,7.65595920e-13,4.21407700e-15
				,-1.01372903e-17,1.18991188e-17,-2.31060000e-18,-1.20159428e-17
				,-1.04083409e-17
		};
		
		for (int i = 0; i < N; i++) {
			System.out.println(fout[i] - p.f[i]);
			assertEquals(fout[i], p.f[i], 5e-10);
		}
/* Python code for generating the above numbers in this test
 
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps

def spectral_to_nodes(c,N):
    c = np.reshape(c,(len(c),1)) #This is to get a column vector, since the fourier transform output depends on the orientation
    cc = np.append(c[1:],np.flip(c[1:-1], axis = 0))/2 #Need to make some changes in order to get the Chebyshev coefficients
    ap = np.append(c[0], cc) #Apply the FFT and throw away the complex values
    F=np.real(np.fft.fft(ap))
    function_values=F[0:N+1] 
    return function_values

def nodes_to_spectral(f,N):
    f1 = np.reshape(f,(len(f),1)) #Need input as a column vector
    f2 = np.reshape(np.flip(f[1:-1], axis = 0),(len(f)-2,1)) #Need to make some alteration before we apply the inverse FFT 
    f = np.real(np.fft.ifft(np.append(f1,f2))) #Apply the inverse fourier transform and throw away the complex values 
    spectral_coefficients = np.append(f[0], 2*f[1:N]) #Get the Chebyshev coefficients
    spectral_coefficients = np.append(spectral_coefficients, f[N])
    return spectral_coefficients

#Example output
n = 20
m = 5
N = 32
xx = N
circle = -np.arange(0,np.pi+np.pi/xx,np.pi/xx)
print("circle=", circle)
x = 1-(np.cos(circle) + 1)/2
print("x=", x)
f = sps.binom.pmf(m, n, x, loc=0)
print("Node values in")
print(f)
plt.plot(x,f)
c = nodes_to_spectral(f,N)
f = spectral_to_nodes(c,N)
print("Chebyshev coefficients")
print(c)
print("Node values out")
print(f)
plt.plot(x,f)
plt.show()
		
 */
	}

	@Test
	public void testChebeshevPolynomial4() {
		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(60, 100);

	// store function values
	double [] f = new double[N];
	System.arraycopy(p.f, 0, f, 0, f.length);
	System.out.println(Arrays.toString(f));
	double [] fin = new double[] {
			0.00000000e+000,9.81217458e-130,8.44740350e-094,5.52746404e-073
			,1.96058626e-058,2.23637075e-047,1.40561393e-038,2.20794099e-031
			,2.12897558e-025,2.24025170e-020,3.79277862e-016,1.35570782e-012
			,1.24328444e-009,3.37020772e-007,2.99255612e-005,9.36366872e-004
			,1.08438667e-002,4.78258544e-002,8.11173785e-002,5.24425537e-002
			,1.25564545e-002,1.05721928e-003,2.88859854e-005,2.27585665e-007
			,4.35710569e-010,1.57796090e-013,7.41815666e-018,2.50114523e-023
			,2.21775979e-030,7.81024821e-040,1.55163732e-053,2.18277324e-077
			,0.00000000e+000
			};
	for (int i = 0; i < N; i++) {
		System.out.println(fin[i] - f[i]);
	}
	double [] a = new double [] {
			6.46372099e-03,-2.55988950e-03,-1.16700845e-02,6.99254795e-03
			,8.30736893e-03,-9.62449434e-03,-3.89177042e-03,9.99997936e-03
			,-3.16413803e-04,-8.41597825e-03,3.32499151e-03,5.70139468e-03
			,-4.70808001e-03,-2.82391503e-03,4.62516343e-03,5.27632266e-04
			,-3.61099170e-03,8.51152558e-04,2.28218230e-03,-1.36076693e-03
			,-1.10635964e-03,1.27694765e-03,3.07749738e-04,-9.21548777e-04
			,1.02340235e-04,5.38834238e-04,-2.32248762e-04,-2.54010441e-04
			,2.18037296e-04,9.22966180e-05,-1.67916650e-04,-2.01820567e-05
			,-1.00910284e-05
	};
	System.out.println(Arrays.toString(p.a));
	for (int i = 0; i < N; i++) {
		System.out.println(a[i] + " " +  p.a[i] + " " + (a[i] - p.a[i]));
//		assertEquals(a[i], p.a[i], 5e-10);
	}
	
	double [] exp = new double[]{
			7.01716030e-03,-2.94255425e-03,-1.11341912e-02,6.33565510e-03
			,5.31269743e-03,-5.42670684e-03,-1.12815346e-03,2.77378802e-03
			,-2.69896780e-04,-9.05679289e-04,2.87824983e-04,1.84696096e-04
			,-1.05203315e-04,-1.94253323e-05,2.28828711e-05,-5.28541986e-07
			,-3.22723820e-06,5.27386625e-07,2.93315541e-07,-8.94030807e-08
			,-1.53099813e-08,8.66637254e-09,1.51203823e-10,-5.45194196e-10
			,4.16486463e-11,2.28086803e-11,-3.75075989e-12,-6.22940046e-13
			,1.69970658e-13,1.01503240e-14,-2.73890476e-15,-7.26617451e-17
			,-6.97734076e-18
	};
	QMatrix Q = new QMatrix(N);
	Q.setQ(new double[]{.5,-1,0,0}, new double[] {0,0.5,-0.5,0});
	double [] b = MatrixExponentiator.expmvRK4(0.1, Q.Q, p.a, 0.001);
	for (int i = 0; i < N; i++) {
		System.out.println((b[i] - exp[i]));
		assertEquals(exp[i], b[i], 5e-8);
	}

	
	
	}
}
