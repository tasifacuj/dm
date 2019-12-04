#include "MainWindow.hpp"

#include <QDebug>

namespace dm {
MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.btn_, &QPushButton::clicked, this, &MainWindow::onClickMeClicked);
}


void MainWindow::onClickMeClicked() {
	qDebug() << __func__;
}
}