package snapper;

import java.util.Arrays;

import snap.likelihood.COMPLEX;
import snap.matrix.AbstractMatrix;

public class QMatrix  extends AbstractMatrix {
	int N;
	double [] D;
	double [] SD;
	double [] S2D;
	double [] S3D;
	double [] D2;
	double [] SD2;
	double [] S2D2;
	double [] Q;

	private double [] x_times_Tx() { //Transpose of Operator B defined as x*f(x) (Gordon's thesis Propositon.5 page. 78)
		    //B = sp.sparse.diags([0.25, 0.5, 0.25], [-1, 0, 1], shape=(N, N)).toarray()
		    //B[1,0] = 0.5
		    double [] B = new double [N*N];
		    for (int i = 0; i < N; i++) {
		    	int j = i - 1;
		    	if (j >= 0) {
		    		B[i * N + j] = 0.25;
		    	}
	    		B[i * N + j + 1] = 0.5;
		    	if (j + 2 < N) {
		    		B[i * N + j + 2] = 0.25;
		    	}
		    }
		    B[N] = 0.5;
		    return B;
	}
	private double [] dot(double [] A, double [] B, double [] S) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double sum = 0;
				for (int k = 0; k < N; k++) {
					sum += A[i*N+k] * B[k*N+j];
				}
				S[i*N+j] = sum;
			}
		}
		return S;		
	}
		private double [] Dd() { // #Derivative operator
		    // N = this.N - 4;
		    //D[0,np.arange(1,N+4,2)] = 2*np.arange(1,N+4,2)
		    //for i in range(0,N+2):
		    //    D[i+1,np.arange(i+2,N+4,2)] = 4*(np.arange(i+1,N+3,2)+1)
		    D = new double [N*N];
		    for (int i = 1; i < N; i++) {
	    		int k = 4 * i;
	    		if (i % 2 == 1) {
	    			D[i] = k / 2;
			    	for (int j = 2; j < i; j += 2) {
			    		D[j*N+i] = k;
			    	}
 	    		} else {
 			    	for (int j = 1; j < i; j += 2) {
 			    		D[j*N+i] = k;
 			    	}
 	    		}
		    }
		    return D;
		}

	
	public QMatrix(int N) {
		this.N = N;
		double [] S = x_times_Tx();
		double [] S2 = new double[N*N];
		dot(S,S, S2);
		double [] S3 = new double[N*N];
		dot(S2,S, S3);
		
		D = Dd();

		
		SD = new double [N*N];
		dot(S,D,SD);
		S2D = new double [N*N];
		dot(S2,D,S2D);
		S3D = new double [N*N];
		dot(S3,D,S3D);

//		prettyPrint(S);
//		prettyPrint(S2);
//		prettyPrint(S3);

		
		D2 = new double [N*N];
		dot(D,D,D2);
		SD2 = new double [N*N];
		dot(S,D2,SD2);
		S2D2 = new double [N*N];
		dot(S2,D2,S2D2);
		
		Q = new double[N*N];
	}

    public void setQ(double [] a, double [] b) {
    	if (Q.length != N * N) {
    		throw new IllegalArgumentException("Exepected matrix of site " + N + "x" + N);
    	}
    	Arrays.fill(Q,  0);
    	add(Q, a[0], D);
    	add(Q, a[1], SD);
    	add(Q, a[2], S2D);
    	add(Q, a[3], S3D);
    	add(Q, 0.5*b[0], D2);
    	add(Q, 0.5*b[1], SD2);
    	add(Q, 0.5*b[2], S2D2);
    }

    public void setQ(double u, double v, double coalescentRate) {
    	if (Q.length != N * N) {
    		throw new IllegalArgumentException("Exepected matrix of site " + N + "x" + N);
    	}
		throw new RuntimeException("Not implemented yet");
    }
    
    
    private void add(double [] Q, double x, double [] M) {
    	if (x == 0) {
    		return;
    	}
		for (int i = 0; i < Q.length; i++) {
			Q[i] += x * M[i];
		}
    }
    
    
    public void prettyPrint(double [] M) {
		System.out.print('[');
    	for (int i = 0; i < N; i++) {
    		if (i > 0) System.out.print(' ');
    		System.out.print('[');
    		for (int j = 0; j < N; j++) {
    			System.out.print(M[i*N+j]);
    			if (j < N-1) {
    				System.out.print(", ");
    			}
    		}
    		if (i == N-1) {
        		System.out.print(']');
    		}
        	System.out.println(']');
    	}
    }
    
	@Override
	public void multiply(double[] v, double[] Av) {
		if (v.length == N + 1) {
			// in case v an Av are indexed [1,..,N] instead of [0,..,N-1]
			for (int i = 0; i < N; i++) {
				double sum = 0;
				for (int k = 0; k < N; k++) {
					sum += Q[i*N+k] * v[k + 1];
				}
				Av[i + 1] = sum;
			}
			return;
		}
		for (int i = 0; i < N; i++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += Q[i*N+k] * v[k];
			}
			Av[i] = sum;
		}
	}
	
	@Override
	public void solve(COMPLEX[] vc, COMPLEX offset, COMPLEX[] xc) throws Exception {
		throw new RuntimeException("not implemented yet");
	}

	@Override
	public void solve(double[] vc_r, double[] vc_i, double offset_r, double offset_i, double[] xc_r, double[] xc_i)
			throws Exception {
		throw new RuntimeException("not implemented yet");		
	}
	
	public void transpose() {
		double [] T = new double[Q.length];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				T[i*N+j] = Q[i+j*N];
			}
		}
		Q = T;		
	}
	
	public int getNrOfRows() {return N;}
	public int getNrOfCols() {return N;}
	public double infNorm() {
		double norm = 0;
		for (double d : Q) {
			norm += Math.abs(d);
		}
		return norm;
	}
	public double trace() {
		double trace = 0;
		for (int i = 0; i < N; i++) {
			trace += Q[i*N+i];
		}
		return trace;
	};
}
