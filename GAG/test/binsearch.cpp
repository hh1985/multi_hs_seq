#include <iostream>
#include <algorithm>

using namespace std;
int main()
{
	double test_data[] = {12.01, 13.57, 15.289, 15.29, 15.292, 18.01, 20.03};
	
	size_t lb = 0, ub = sizeof(test_data)/sizeof(double), mid;
	double ppm1 = -0.01;
	double ppm2 = 0.01;

	double ms = 15.295; 

	for(; lb < ub; )
	{
		mid = (lb+ub)/2;
		double distance = ms - test_data[mid];
		if(distance > ppm2) {
			ub = mid - 1;
		} else if(distance < ppm1) {
			lb = mid + 1;
		} else {
			cout << "Found " << test_data[mid] << endl;
			// Expand.
			// size_t cur = mid;
			for(size_t cur = mid+1; cur <= ub; cur++)
			{
				distance = ms - test_data[cur];
				if(distance <= ppm2 && distance >= ppm1)
					cout << "Found " << test_data[cur] << endl;
				else
					break;
			}
			for(size_t cur = mid-1; cur >= lb; cur--)
			{
				distance = ms - test_data[cur];
				if(distance <= ppm2 && distance >= ppm1)
					cout << "Found " << test_data[cur] << endl;
				else
					break;
			}
			break;
		}
	}
}