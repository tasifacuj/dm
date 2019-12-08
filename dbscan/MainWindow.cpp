#include "MainWindow.hpp"
#include "kd_tree/PointCloud.hpp"
#include "kd_tree/kd-tree.hpp"
#include "dbscan/PointCloud2D.hpp"
#include "dbscan/DbScan.hpp"

#include <QDebug>
#include <QFile>
#include <QTextStream>

#include <cstdlib>
#include <array>


namespace dm {



template <typename T>
void generateRandomPointCloud(KdTree::PointCloud<T> &point, const size_t N, const T max_range = 10)
{
	// Generating Random Point Cloud
	point.pts.resize(N);
	
	for (size_t i = 0; i < N; i++){
		point.pts[i].x = max_range * (rand() % 1000) / T(1000);
		point.pts[i].y = max_range * (rand() % 1000) / T(1000);
		point.pts[i].z = max_range * (rand() % 1000) / T(1000);
	}
}

template <typename T>
void generateRandomPointCloud2D(dbscan::PointCloud2D<T> &points, const size_t N, const T max_range = 10, float offset = 0.0f){
	for (size_t i = 0; i < N; i++) {
		typename dbscan::PointCloud2D<T>::Point point{ max_range * (rand() % 1000) / T(1000) + offset, max_range * (rand() % 1000) / T(1000) + offset, dbscan::PointCloud2D<T>::Undefined };
		points.pts.push_back(point);
	}
}

template <typename num_t>
void kdtree_demo(const size_t N){
	KdTree::PointCloud<num_t> cloud;

	// Generate points:
	generateRandomPointCloud(cloud, N);
	//cloud.pts.push_back(dbscan::PointCloud2D<num_t>::Point{ 13.5f, 13.5f, -1 });
	// construct a kd-tree index:
	typedef KdTree::KDTreeSingleIndexAdaptor<
		KdTree::L2_Simple_Adaptor<num_t, KdTree::PointCloud<num_t> >,
		KdTree::PointCloud<num_t>,
		3 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KdTree::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();
	const num_t query_pt[3] = { cloud.pts[ 1].x, cloud.pts[1].y, cloud.pts[1].z };
	// ----------------------------------------------------------------
	// radiusSearch(): Perform a search for the points within search_radius
	// ----------------------------------------------------------------
	{
		const num_t search_radius = static_cast<num_t>(0.1);
		std::map< size_t, num_t >   ret_matches;

		KdTree::SearchParams params;
		//params.sorted = false;

		const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

		qDebug() << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches";
		assert(nMatches == ret_matches.size());
		auto it = ret_matches.begin();

		for (size_t i = 0 ; it != ret_matches.end(); it++, i++)
			qDebug() << "idx[" << i << "]=" << it->first << " dist[" << i << "]=" << it->second;
	}
}

template <typename num_t>
void MainWindow::dbscanDemo(const size_t N, QChartView& chartView) {
	dbscan::PointCloud2D<num_t> cloud;

	// Generate points:
	generateRandomPointCloud2D(cloud, N, 1.5f, 0.0f);
	generateRandomPointCloud2D(cloud, N, 1.5f, 5.0f);
	generateRandomPointCloud2D(cloud, 10, 10.0f, 0.0f);
	
	typedef KdTree::KDTreeSingleIndexAdaptor<
		KdTree::L2_Simple_Adaptor<num_t, dbscan::PointCloud2D<num_t> >,
		dbscan::PointCloud2D<num_t>,
		2 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(2 /*dim*/, cloud, KdTree::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();
	
	// ----------------------------------------------------------------
	// dbscan
	// ----------------------------------------------------------------
	{
		typedef dbscan::DbScan< dbscan::PointCloud2D<num_t>, my_kd_tree_t, size_t> DbScanT;
		DbScanT dbscan(cloud, index);
		dbscan::ScanParams scanParams(0.2, 5);
		dbscan.evaluate( scanParams );

		for (size_t idx = 0; idx < cloud.size(); idx++) {
			qDebug() << "p[" << idx << "] = {" << cloud.pts[idx].x << ", " << cloud.pts[idx].y << "}, cluster id = " << cloud.pts[idx].label;
		}	

		std::array<QColor, 10> colors{
			QColor::fromRgb(0,0,0),
			QColor::fromRgb(34, 193, 40),
			QColor::fromRgb(61, 191, 210),
			QColor::fromRgb(217, 223, 30),
			QColor::fromRgb(255,51,51),
			QColor::fromRgb(255,102,255),
			QColor::fromRgb(255,0,0),
			QColor::fromRgb( 169, 179, 74 ),
			QColor::fromRgb( 74, 119, 179 ),
			QColor::fromRgb( 167, 188, 62 )
		};

		auto clusterIds = cloud.getClusterIds();
		int i = 0;
		for ( auto it = clusterIds.begin(), it_end = clusterIds.end(); it != it_end; it++, i++ ){
			QScatterSeries* series = nullptr;
			if ( series_.count(*it)) {
				auto s_it = series_.find(*it);
				assert(s_it != series_.end());
				series = s_it->second;
				series->clear();
			}else {
				series = new QScatterSeries;
				series->setName(QString("cluster %1").arg(*it));
				series->setMarkerShape(QScatterSeries::MarkerShapeCircle);
				series->setMarkerSize(10);
				series->setBrush(colors[i]);
				series_.emplace( *it, series);
			}
			auto clusterPoints = cloud.getPointsByClusterId(*it);
			
			for (auto p : clusterPoints)
				*series << QPointF(p.x, p.y);

			chartView.chart()->addSeries(series);
			chartView.chart()->setTitle("DbScan example");
			chartView.chart()->createDefaultAxes();
			chartView.chart()->setDropShadowEnabled(false);
			chartView.chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
			QAbstractAxis* x = chartView.chart()->axisX();
			x->setRange(-1, 10);
			QAbstractAxis* y = chartView.chart()->axisY();
			y->setRange(-1, 10);
		}
	}
}


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.btn_, &QPushButton::clicked, this, &MainWindow::onClickMeClicked);
	QObject::connect(ui.loadFileBtn_, &QPushButton::clicked, this, &MainWindow::onLoadClicked);
	srand(static_cast<unsigned int>(time(nullptr)));
}


void MainWindow::onClickMeClicked() {

	//kdtree_demo<float>(10);
	//kdtree_demo<double>(100000);
	dbscanDemo<float>(50, *ui.graphicsView);
	qDebug() << __func__;
	
}

void MainWindow::onLoadClicked(){
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open Points"));
	qDebug() << "Selected " << fileName;
	QFile file(fileName);
	
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		qCritical() << "Can't open " << fileName;

	QTextStream in(&file);
	dbscan::PointCloud2D<float> cloud;

	while (!in.atEnd()) {
		QString line = in.readLine(); //read one line at a time
		QStringList cols = line.split(" ");

		if (cols.size() != 2)
			continue;

		dbscan::PointCloud2D<float>::Point point{ cols[ 0 ].toDouble(), cols[ 1 ].toDouble(), dbscan::PointCloud2D<float>::Undefined };
		cloud.pts.push_back(point);
	}
	qDebug() << "Got " << cloud.pts.size() << " points";
	scan<dbscan::PointCloud2D<float>,float>(cloud, *ui.graphicsView);
}

template <typename Cloud, typename num_t>
void MainWindow::scan(Cloud& cloud, QChartView& chartView) {
	typedef KdTree::KDTreeSingleIndexAdaptor<
		KdTree::L2_Simple_Adaptor<num_t, dbscan::PointCloud2D<num_t> >,
		dbscan::PointCloud2D<num_t>,
		2 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(2 /*dim*/, cloud, KdTree::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();

	// ----------------------------------------------------------------
	// dbscan
	// ----------------------------------------------------------------
	{
		typedef dbscan::DbScan< dbscan::PointCloud2D<num_t>, my_kd_tree_t, size_t> DbScanT;
		DbScanT dbscan(cloud, index);
		dbscan::ScanParams scanParams(10, 15);
		dbscan.evaluate(scanParams);
		int clusters = cloud.getClusterNum();
		qDebug() << "Got " << clusters << " clusters";
		for (size_t idx = 0; idx < cloud.size(); idx++) {
			qDebug() << "p[" << idx << "] = {" << cloud.pts[idx].x << ", " << cloud.pts[idx].y << "}, cluster id = " << cloud.pts[idx].label;
		}

		std::array<QColor, 20> colors{
			QColor::fromRgb(0,0,0),
			QColor::fromRgb(34, 193, 40),
			QColor::fromRgb(61, 191, 210),
			QColor::fromRgb(217, 223, 30),
			QColor::fromRgb(255,51,51),
			QColor::fromRgb(255,102,255),
			QColor::fromRgb(255,0,0),
			QColor::fromRgb(169, 179, 74),
			QColor::fromRgb(74, 119, 179),
			QColor::fromRgb(167, 188, 62),
			QColor::fromRgb(0,0,0),
			QColor::fromRgb(34, 193, 40),
			QColor::fromRgb(61, 191, 210),
			QColor::fromRgb(217, 223, 30),
			QColor::fromRgb(255,51,51),
			QColor::fromRgb(255,102,255),
			QColor::fromRgb(255,0,0),
			QColor::fromRgb(169, 179, 74),
			QColor::fromRgb(74, 119, 179),
			QColor::fromRgb(167, 188, 62)
		};

		auto clusterIds = cloud.getClusterIds();
		int i = 0;
		for (auto it = clusterIds.begin(), it_end = clusterIds.end(); it != it_end; it++, i++) {
			QScatterSeries* series = nullptr;

			if (i == 20) {
				qWarning() << "10 clusters reached";
				return;
			}

			if (series_.count(*it)) {
				auto s_it = series_.find(*it);
				assert(s_it != series_.end());
				series = s_it->second;
				series->clear();
			}
			else {
				series = new QScatterSeries;
				series->setName(QString("cluster %1").arg(*it));
				series->setMarkerShape(QScatterSeries::MarkerShapeCircle);
				series->setMarkerSize(10);
				series->setBrush(colors[i]);
				series_.emplace(*it, series);
			}
			auto clusterPoints = cloud.getPointsByClusterId(*it);

			for (auto p : clusterPoints)
				*series << QPointF(p.x, p.y);

			chartView.chart()->addSeries(series);
			chartView.chart()->setTitle("DbScan example");
			chartView.chart()->createDefaultAxes();
			chartView.chart()->setDropShadowEnabled(false);
			chartView.chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
			QAbstractAxis* x = chartView.chart()->axisX();
			x->setRange(-1, 500);
			QAbstractAxis* y = chartView.chart()->axisY();
			y->setRange(-1, 700);
		}
	}
}

}// namespace dm