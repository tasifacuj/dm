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
		void onLoadClicked();
	private:
		template <typename num_t>
		void dbscanDemo(const size_t N, QChartView& chartView);

		template <typename Cloud, typename num_t>
		void scan(Cloud& cloud, QChartView& chartView);
	private:
		Ui::MainWindowClass ui;
		std::map< int, QScatterSeries* > series_;
	};
}
