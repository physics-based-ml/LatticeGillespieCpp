#include <iostream>
#include <map>
#include <numeric>
#include <vector>
#include <ctime>

struct Site {
	int x;
	int y;
	int z;
};
bool operator <(const Site& a, const Site& b) {
	return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
};

struct Mol {
	std::string s;
};

int main() {

	std::map<Site,Mol> m;
	std::map<int,std::map<int,std::map<int,Mol>>> n;

	// Vectors of sites
	std::vector<int> x1(30);
	std::iota(std::begin(x1), std::end(x1), 1);
	std::vector<int> x2 = x1, x3 = x1;

	// Shuffle
	std::random_shuffle(x1.begin(),x1.end());
	std::random_shuffle(x2.begin(),x2.end());
	std::random_shuffle(x3.begin(),x3.end());

	// Find a random free site
	std::map<Site,Mol>::iterator itmap;
	Site s;
	Mol mol;
	int ctr=0;
	for (auto i1=0; i1 < 15; i1++) {
		for (auto i2=0; i2 < 15; i2++) {
			for (auto i3=0; i3 < 15; i3++) {
				s.x = x1[i1];
				s.y = x2[i2];
				s.z = x3[i3];
				mol.s = "o";
				m.insert(std::make_pair(s,mol));
				n[x1[i1]][x2[i2]][x3[i3]] = mol;
			};
		};
	};

	// Search
	std::clock_t    start;

	start = std::clock();
	s.x = x1[10];
	s.y = x2[10];
	s.z = x3[10];
	auto it = m.find(s);
	if (it != m.end()) {
		std::cout << it->second.s << std::endl;
	};
	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	start = std::clock();
	s.x = 80;
	s.y = x2[10];
	s.z = x3[10];
	it = m.find(s);
	if (it != m.end()) {
		std::cout << it->second.s << std::endl;
	};
	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	start = std::clock();
	std::map<int,std::map<int,std::map<int,Mol>>>::iterator it1;
	std::map<int,std::map<int,Mol>>::iterator it2;
	std::map<int,Mol>::iterator it3;
	it1 = n.find(x1[10]);
	if (it1 != n.end()) {
		it2 = it1->second.find(x2[10]);
		if (it2 != it1->second.end()) {
			it3 = it2->second.find(x3[10]);
			if (it3 != it2->second.end()) {
				std::cout << it3->second.s << std::endl;
			};
		};
	};
	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	start = std::clock();
	it1 = n.find(80);
	if (it1 != n.end()) {
		it2 = it1->second.find(x2[10]);
		if (it2 != it1->second.end()) {
			it3 = it2->second.find(x3[10]);
			if (it3 != it2->second.end()) {
				std::cout << it3->second.s << std::endl;
			};
		};
	};
	std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

	return 0;
};