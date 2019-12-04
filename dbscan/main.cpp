#include "MainWindow.hpp"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	dm::MainWindow w;
	w.show();
	return a.exec();
}
