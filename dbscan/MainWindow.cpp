#include "MainWindow.hpp"
#include "kd_tree/PointCloud.hpp"
#include "kd_tree/kd-tree.hpp"
#include "dbscan/PointCloud2D.hpp"
#include "dbscan/DbScan.hpp"

#include <QDebug>

#include <cstdlib>

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
void dbscan_demo(const size_t N) {
	dbscan::PointCloud2D<num_t> cloud;

	// Generate points:
	generateRandomPointCloud2D(cloud, N, 1.0f, 0.0f);
	generateRandomPointCloud2D(cloud, N, 1.0f, 10.0f);
	
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
		dbscan::ScanParams scanParams(0.2, 10);
		dbscan.evaluate( scanParams );

		for (size_t idx = 0; idx < cloud.size(); idx++) {
			qDebug() << "p[" << idx << "] = {" << cloud.pts[idx].x << ", " << cloud.pts[idx].y << "}, cluster id = " << cloud.pts[idx].label;
		}
	}
}


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	QObject::connect(ui.btn_, &QPushButton::clicked, this, &MainWindow::onClickMeClicked);
	srand(static_cast<unsigned int>(time(nullptr)));
}


void MainWindow::onClickMeClicked() {

	//kdtree_demo<float>(10);
	//kdtree_demo<double>(100000);
	dbscan_demo<float>(50);
	qDebug() << __func__;
}

}