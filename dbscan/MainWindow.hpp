#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_MainWindow.h"

namespace dm {
	class MainWindow : public QMainWindow
	{
		Q_OBJECT

	public:
		MainWindow(QWidget *parent = Q_NULLPTR);
	private slots:
		void onClickMeClicked();
	private:
		Ui::MainWindowClass ui;
	};
}
