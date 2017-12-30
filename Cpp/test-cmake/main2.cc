#include <boost/smart_ptr/shared_ptr.hpp>
#include <iostream>
#include <memory>

class A{};

class classA
{
	boost::shared_ptr<int> ptA;
public:
	classA(boost::shared_ptr<int> p) : ptA(p) {}
	void setValue(int n) {
		*ptA = n;
	}
};

class classB
{
	boost::shared_ptr<int> ptB;
public:
	classB(boost::shared_ptr<int> p) : ptB(p) {}
	int getValue() const {
		return *ptB;
	}
};

int main()
{
	boost::shared_ptr<A> pA(new A);

	std::cout << pA.get() << std::endl;

	boost::shared_ptr<A> pB(pA);

	std::cout << pA.get() << std::endl;
	std::cout << pB.get() << std::endl;
	
	
	boost::shared_ptr<int> pTemp(new int(2013));
	classA a(pTemp);
	classB b(pTemp);

	a.setValue(2014);
	std::cout << "b.getValue() = " << b.getValue() << std::endl;

	return 0;
}