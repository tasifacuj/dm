#pragma once

#include <QtWidgets/QMainWindow>
#include <QScatterSeries>
#include "ui_MainWindow.h"

#include <set>


namespace dm {
	class MainWindow : public QMainWindow
	{
		Q_OBJECT

	public:
		MainWindow(QWidget *parent = Q_NULLPTR);
	private slots:
		void onClickMeClicked();
	private:
		template <typename num_t>
		void dbscanDemo(const size_t N, QChartView& chartView);
	private:
		Ui::MainWindowClass ui;
		std::map< int, QScatterSeries* > series_;
	};
}
