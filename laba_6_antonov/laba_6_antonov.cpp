
#include <iostream>
#include <mpi.h>
#include <fstream>
//#include <math.h>
#include <cmath>
#include <stdlib.h>


const int root = 0, tag = 0;
const int n = 4000, m = 4000;
double A[n][m], B[m][n], v[n], d[n], C[n][n], A1[n][m], B1[m][n], v1[n];

void FillMatrix(double(&AA)[n][m], double(&BB)[m][n]) {
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			AA[i][j] = 1/*random() % 100 / 2.4*/;
		}
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			BB[i][j] = 1/*random() % 100 / 2.4*/;
		}
	}
}
/*
void FillVector(double v1[n]) {
	for (int j = 0; j < n; j++) {
		v1[j] = 1/*rand() % 100 / 2.4;
	}
}
*/

void Matrix_Peremnoj(double(&AA)[n][m], double(&BB)[m][n], double(&CC)[n][n]) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			CC[i][j] = 0;
			for (int k = 0; k < m; k++) {
				CC[i][j] += AA[i][k] * BB[k][j];
			}
		}
	}
}
void Matrix_Peremnoj_na_vector(double(&AA)[n][m], double(&vv)[n]) {
	for (int i = 0; i < n; i++) {
		d[i] = 0;
		for (int j = 0; j < m; j++) {
			d[i] += vv[i] * AA[i][j];
		}
	}
}


void Zapis_v_File() {
	std::ofstream File1("./Matrix_1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			File1 << A[i][j] << " ";
		}
		File1 << std:: endl;
	}
	File1.close();

	std::ofstream File2("./Matrix_2.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			File2 << B[i][j] << " ";
		}
		File2 << std::endl;
	}
	File2.close();

	//ofstream File4("C:/Users/neste/source/Repos/laba_6_antonov/Vector_1.txt");
	//for (int i = 0; i < n; i++)
	//{
	//	File4 << v[i] << endl;
	//}
	//File4.close();
}

void Zapix_otvetov_v_File(double(&CC)[n][n]/*,double(&d)[n]*/) {
	std::cout << " c zapisannaya --" << std::endl;
	std::ofstream File3("./Matrix_Otvet1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			File3 << CC[i][j] << " ";
		}
		File3 << "\n";
	}
	File3.close();
	
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++) {
	//		std::cout << CC[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	/*ofstream File5("/home/vc/18VF1/laba2-master/laba2_try3/Vector_Otvet1.txt");
	for (int i = 0; i < n; i++)
	{
		File5 << d[i] << endl;
	}
	File5.close();*/
}

//void read_Vector() {
//	ifstream File5("C:/Users/neste/source/Repos/laba_6_antonov/Vector_1.txt");
//	for (int i = 0; i < n; i++) {
//		File5 >> v1[i];
//		//  cout << DD[i]<<endl;
//	}
//	File5.close();
//}
void read_MatrixB() {
	std::ifstream File2("./Matrix_2.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			File2 >> B1[i][j];
			//cout << AA[i][j] << " ";
		}
		// cout << endl;
	}
	File2.close();
}
void read_MatrixA() {
	std::ifstream File1("./Matrix_1.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++) {
			File1 >> A1[i][j];
			//cout << AA[i][j] << " ";
		}
		// cout << endl;
	}
	File1.close();

}

double k[m],K[n][m];
double l[m],L[n][m];

void vzat_vector_iz_matrixA(double(&AA)[n][m], int s1/* double(&BB)[m][n], int s2,*/ ) {
	// if c == 1 znachit berem stroku if c== 2 znachit berem stolbec
	// s eto nomer stroki ili stolbca kotoriy nuzhno vzat

	for (int j = 0; j < m; j++) {
		k[j] = AA[s1][j];
		//cout << k[j] << endl;
	}
}
void vzat_minimatrix_iz_matrixA(double(&AA)[n][m], int s1, int s2/* double(&BB)[m][n], int s2,*/) {
	// if c == 1 znachit berem stroku if c== 2 znachit berem stolbec
	// s eto nomer stroki ili stolbca kotoriy nuzhno vzat
	for (int i = s1; i < s2; i++)
	{
		for (int j = 0; j < m; j++) {
			K[i][j] = AA[i][j];
			//cout << k[j] << endl;
		}
	}
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < m; j++) {
	//		std::cout << K[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}
}
void vzat_minimatrix_iz_matrixB(double(&BB)[n][m], int s1, int s2/* double(&BB)[m][n], int s2,*/) {
	// if c == 1 znachit berem stroku if c== 2 znachit berem stolbec
	// s eto nomer stroki ili stolbca kotoriy nuzhno vzat
	for (int i = s1; i < s2; i++)
	{
		for (int j = 0; j < m; j++) {
			L[j][i] = BB[j][s1];
			//cout << k[j] << endl;
		}
	}
}
void vzat_vector_iz_matrixB(double(&BB)[m][n], int s1/* double(&BB)[m][n], int s2,*/) {
	// if c == 1 znachit berem stroku if c== 2 znachit berem stolbec
	// s eto nomer stroki ili stolbca kotoriy nuzhno vzat

	for (int j = 0; j < n; j++) {
		l[j] = BB[j][s1];
		//cout << l[j] << endl;
	}
	
}

double peremnoj_vector_na_vector(double(&kk)[m], double(&ll)[m]) {
	double h = 0;
	for (int i = 0; i < m; i++) {
		h += kk[i] * ll[i];
	}
	return h;
}

double Temp[n][n];
int Mesto[] = { 0,0 };
int N = 0, M = 0, Nn = 0, Mm = 0, limit_1 = 0;
//int kol_strok_v_posled_str_bloke = 0;
//int kol_stolb_v_posled_stolb_bloke = 0;
//bool is_n_menwe_size;




double  rbufA[n][m], rbufB[m][n], /*int gsize; int asize; int bsize;*/  buff[n][n]{};
int main() {
	MPI_Init(NULL, NULL);
	//std::cout << " ya der'mo1" << std::endl;
	double starttime, endtime;
	starttime = MPI_Wtime();

	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int color = world_rank / world_size;

	int popa = n / world_size;
	int pisa = n % world_size;
	

	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

	int row_rank, row_size;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_size(row_comm, &row_size);

	//
	//int nmin = 0, nextra = 0;
	//
	//std::cout << " ya der'mo" << std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);


	//

	double readtime0 = MPI_Wtime();
	if (row_rank == 0) {
		//FillMatrix(A, B);
		//FillVector(v);
		//Zapis_v_File();
		//read_Vector();
		//read_Matrix();
		
		read_MatrixA();
	}
	if (row_rank == 1) {
		read_MatrixB();
	}
	//double readtimer = MPI_Wtime();
	MPI_Bcast(A1, n * m, MPI_DOUBLE, 0, row_comm);
	MPI_Bcast(B1, n * m, MPI_DOUBLE, 1, row_comm);
	//double readtimer1 = MPI_Wtime();
	//std::cout << "read time A= " << readtimer1 - readtimer << std::endl;

	double readtime1 = MPI_Wtime();
	std::cout << "read time B= " << readtime1 - readtime0 << std::endl;
	/*double readtime2 = MPI_Wtime();
	MPI_Bcast(A1, n * m, MPI_DOUBLE, 0, row_comm);
	double readtime3 = MPI_Wtime();
	std::cout << "bcast time = " << readtime3 - readtime2 << std::endl;*/
	
	////int sendcounts[n], displa[n];
	////MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	//
	////MPI_Bcast(B1, n * m, MPI_DOUBLE, 0, row_comm);
	////MPI_Bcast(A1, n * m, MPI_DOUBLE, 0, row_comm);

	


	for (int i = 0; i < n - (n % row_size); i += n / row_size)
	{
		//vzat_vector_iz_matrixA(A1, i);
		//if (i % 2 == 0){
		//cout << "jopa" <<endl;
		
		/*if (row_rank == 0)
		{*/
			//vzat_minimatrix_iz_matrixA(A1, i, i + n/row_size);
		//}
		/*double readtime2 = MPI_Wtime();
		MPI_Bcast(K, n/row_size * m, MPI_DOUBLE, 0, row_comm);
		double readtime3 = MPI_Wtime();
		std::cout <<"rank = "<<row_rank<< " bcast time = " << readtime3 - readtime2 << std::endl;*/

		//MPI_Scatter(K,  row_size, MPI_DOUBLE, rbufA,  row_size, MPI_DOUBLE, 0, row_comm);
		//MPI_Scatterv(k, sendcounts, displa, MPI_DOUBLE, rbufA, n / row_size, MPI_DOUBLE, 0, row_comm);
		//MPI_Scatter(A1, n / row_size, MPI_DOUBLE, rbufA, n / row_size, MPI_DOUBLE, 0, row_comm);

		//}else{
		//}
		for (int j = 0; j < m; j += m / row_size)
		{
			//vzat_vector_iz_matrixB(B1, j);
			//if (i % 2 == 0){
			//vzat_minimatrix_iz_matrixB(B1, j, j + row_size - 1);
			//cout << sendcounts[i] << " " << displa[i] << endl;
			//MPI_Scatterv(l, sendcounts, displa, MPI_DOUBLE, rbufB, n / row_size, MPI_DOUBLE, 0, row_comm);
			//MPI_Scatter(B1, n * m  , MPI_DOUBLE, rbufB, n * m , MPI_DOUBLE, 1, row_comm);
			double readtimer = MPI_Wtime();

			for (int ii = row_rank * (n / row_size); ii < (row_rank + 1) * (n / row_size); ii++) {
				for (int jj = row_rank * (n / row_size); jj < (row_rank + 1) * (n / row_size); jj++) {
					//CC[i][j] = 0;
					for (int kk = 0; kk < m; kk++) {
						Temp[i + ii][j + jj] += A1[i+ii][kk] * B1[kk][j+jj];

					}//if (row_rank == 1)std::cout << Temp[i + ii][j + jj];
				}//if (row_rank == 1)std::cout<<std::endl;
			}
			double readtimer1 = MPI_Wtime();
			std::cout << "count time = " << readtimer1 - readtimer <<" rank = "<< row_rank<<" quad #"<< i+j<< std::endl;
			//}else {
			//MPI_Scatter(l, 1, MPI_DOUBLE, rbufB, 1, MPI_DOUBLE, 1, comm_b);
			//}
			//Matrix_Peremnoj(K, rbufB, Temp);
			//for (int ii = row_rank * n / row_size; ii < (row_rank + 1) * n / row_size; ii++) {
			//	for (int jj = row_rank * n / row_size; jj < (row_rank + 1) * n / row_size; jj++) {
			//		//CC[i][j] = 0;
			//		for (int kk = 0; kk < m; kk++) {
			//			Temp[i + ii][j + jj] += A1[i+ii][j+kk] * rbufB[kk][jj];
			//		}
			//	}
			//}
			/*for (int i = 0; i < n / row_size; i++)
			{
				Temp[i] = rbufA[i] * rbufB[i];
			}*/
			//Temp[0] = rbufA[0] * rbufB[0];
			//cout << rbufA[0] << " " << rbufB[0] << " " << rank << endl;
			//if (i % 2 == 0){
			//MPI_Reduce(Temp, buff, n / row_size, MPI_DOUBLE, MPI_SUM, 0, row_comm);
			// 
			MPI_Gather(&Temp[i][j], n * m / row_size/row_size, MPI_DOUBLE, &buff[i][j], n * m / row_size/row_size, MPI_DOUBLE, 0, row_comm);
			// 
			//}else{
			//MPI_Reduce(Temp, buff, 1, MPI_DOUBLE, MPI_SUM, 1, comm_b);
			//}
			//cout << buff[0] << " " << rank << endl;
			//if (row_rank == 0)
			//{
			//	//cout << buff[0] << " " << rank << endl;
			//	for (int i = 1; i < n/row_size; i++)
			//	{
			//		buff[0] += buff[i];
			//	}
			//	C[i][j] = buff[0];
			//}

			fflush(stdout);
		}
	}
	//if (n % row_size > 0) {
	//	for (int i = 0; i < n - (n % row_size); i += row_size)
	//	{
	//		vzat_minimatrix_iz_matrixA(A1, i, i + row_size - 1);
	//		for (int j = 0; j < m; j += world_size)
	//		{
	//			vzat_minimatrix_iz_matrixB(B1, j, j + row_size - 1);
	//			MPI_Scatter(L, row_size, MPI_DOUBLE, rbufB, row_size, MPI_DOUBLE, 0, row_comm);
	//			Matrix_Peremnoj(K, rbufB, Temp);
	//			MPI_Gather(Temp, row_size, MPI_DOUBLE, &buff[i][j], row_size, MPI_DOUBLE, 0, row_comm);
	//			fflush(stdout);
	//		}
	//	}
	//}
	MPI_Barrier(row_comm);
	if (row_rank == 0) {

		
		double starttimeZ = MPI_Wtime();
		Zapix_otvetov_v_File(buff/* buff, d */ );
		endtime = MPI_Wtime();
		printf("vipolnenie zanyalo %f seconds\n", endtime - starttime);
		printf("Zapis zanyala %f seconds\n", endtime - starttimeZ);
	}
	MPI_Barrier(row_comm);

	MPI_Comm_free(&row_comm);
	MPI_Finalize();
	////MPI_Finalize();
	return 1;
}
