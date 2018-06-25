
#include "sim_sop.h"
#include <SYS/SYS_Math.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_Node.h>
#include <UT/UT_DSOVersion.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <PRM/PRM_Template.h>


static PRM_Name X("x", "x");
static PRM_Name Y("y", "y");
static PRM_Name timeG("time", "time");
static PRM_Name iteration("iteration", "iteration");


OP_Node * Simple_SOP::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {

	return new Simple_SOP(net, name, op);


		

}
PRM_Template Simple_SOP::MyTemplateList[] = {


	PRM_Template(PRM_XYZ,1, &X),
	PRM_Template(PRM_XYZ,1, &Y),
	PRM_Template(PRM_XYZ,1, &timeG),
	PRM_Template(PRM_XYZ,1, &iteration),
	PRM_Template(),


};
Simple_SOP::Simple_SOP(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op) {}
Simple_SOP::~Simple_SOP() {}





void updateWave(fpreal timeInterval, fpreal* x, fpreal* speed) {
	(*x) += timeInterval * (*speed);

	// Boundary reflection
	if ((*x) > 1.0) {
		(*speed) *= -1.0;
		(*x) = 1.0 + timeInterval * (*speed);

	}
	else if ((*x) < 0.0) {
		(*speed) *= -1.0;
		(*x) = timeInterval * (*speed);

	}
}


OP_ERROR
Simple_SOP::cookMySop(OP_Context &context) {

	OP_Node::flags().timeDep = 1;
	fpreal test = context.getTime();
	test *= 0.1;

	GA_Attribute *zestx = gdp->addFloatTuple(GA_ATTRIB_DETAIL, "x", 1);
	GA_Attribute *zesty = gdp->addFloatTuple(GA_ATTRIB_DETAIL, "y", 1);

	GA_RWHandleF h(zestx);
	GA_RWHandleF h1(zesty);
	GA_Offset offset = 0;



	double x = 0;
	double y = 0;
	fpreal speedx = xx(test);
	fpreal speedy = yy(test);

	speedx /= 100 * TIMES(test);
	speedy /= 100 * TIMES(test);


	for (int i = 0; i < Iterations(test); ++i) {
		updateWave(test, &x, &speedx);
		updateWave(test, &y, &speedy);
	}



	h.set(0, fpreal(x));
	h1.set(0, fpreal(y));

	return error();
}








void
newSopOperator(OP_OperatorTable *table)
{
	OP_Operator *op;
	op = new OP_Operator("timetest",
		"timetest",
		Simple_SOP::myConstructor,
		Simple_SOP::MyTemplateList,
		0,
		0,
		0,
		OP_FLAG_GENERATOR);

	table->addOperator(op);
}



