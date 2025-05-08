// Copyright (C) 2004 by Olivier Galizzi, Janek Kozicki                  *
// © 2008 Václav Šmilauer
#pragma once

#ifndef Q_MOC_RUN
#include <core/Omega.hpp>
#include <pkg/common/OpenGLRenderer.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#endif

#include <QGLViewer/constraint.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>
#include <QOpenGLWidget>

namespace yade { // Cannot have #include directive inside.

/*! Class handling user interaction with the openGL rendering of simulation.
 *
 * Clipping planes:
 * ================
 *
 * Clipping plane is manipulated after hitting F1, F2, .... To end the manipulation, press Escape.
 *
 * Keystrokes during clipping plane manipulation:
 * * space activates/deactives the clipping plane
 * * x,y,z aligns the plane with yz, xz, xy planes
 * * left-double-click aligns the plane with world coordinates system (canonical planes + 45˚ interpositions)
 * * 1,2,... will align the current plane with #1, #2, ... (same orientation)
 * * r reverses the plane (normal*=-1)a
 *
 * Keystrokes that work regardless of whether a clipping plane is being manipulated:
 * * Alt-1,Alt-2,... adds/removes the respective plane to bound group:
 * 	mutual positions+orientations of planes in the group are maintained when one of those planes is manipulated
 *
 * Clip plane number is 3; change YADE_RENDERER_NUM_CLIP_PLANE, complete switches "|| ..." in keyPressEvent
 * and recompile to have more.
 */
class GLViewer : public QGLViewer {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wsuggest-override"
	Q_OBJECT
#pragma GCC diagnostic pop

	friend class QGLThread;

protected:
	shared_ptr<OpenGLRenderer> renderer;

private:
	bool                                   isMoving;
	bool                                   wasDynamic;
	float                                  cut_plane;
	int                                    cut_plane_delta;
	bool                                   gridSubdivide;
	bool                                   displayGridNumbers;
	bool                                   autoGrid;
	Real                                   prevGridStep;
	Real                                   requestedGridStep;
	int                                    prevSegments;
	int                                    gridDecimalPlaces;
	Vector3r                               gridOrigin;
	long                                   last;
	int                                    manipulatedClipPlane;
	std::set<int>                          boundClipPlanes;
	shared_ptr<qglviewer::LocalConstraint> xyPlaneConstraint;
	string                                 strBoundGroup()
	{
		string ret;
		FOREACH(int i, boundClipPlanes) ret += " " + boost::lexical_cast<string>(i + 1);
		return ret;
	}
	boost::posix_time::ptime last_user_event;

	static void staticCloseEvent(QCloseEvent* e, const int);

public:
	void updateGLViewer()
	{
#if QGLVIEWER_VERSION < 0x020700
		this->updateGL();
#else
		this->update();
#endif
	};

	const int viewId;

	void centerMedianQuartile();
	int  drawGrid;
	bool drawScale;
	int  timeDispMask;
	enum { TIME_REAL = 1, TIME_VIRT = 2, TIME_ITER = 4 };

	GLViewer(int viewId, const shared_ptr<OpenGLRenderer>& renderer, QOpenGLWidget* = 0);
	virtual ~GLViewer();
#if 0
			virtual void paintGL();
#endif
	void draw() override;
	void drawWithNames() override;
	void displayMessage(const std::string& s) { QGLViewer::displayMessage(QString(s.c_str())); }
	void centerScene(
	        const Real&     suggestedRadius      = -1,
	        const Vector3r& setGridOrigin        = Vector3r(0, 0, 0),
	        const Vector3r& suggestedCenter      = Vector3r(0, 0, 0),
	        int             setGridDecimalPlaces = 4);
	void     centerPeriodic();
	void     mouseMovesCamera();
	void     mouseMovesManipulatedFrame(qglviewer::Constraint* c = NULL);
	void     resetManipulation();
	bool     isManipulating();
	void     startClipPlaneManipulation(int planeNo);
	Real     getSuggestedRadius() const;
	Vector3r getGridOrigin() const;
	Vector3r getSuggestedCenter() const;
	int      getGridDecimalPlaces() const;
	//! get QGLViewer state as string (XML); QGLViewer normally only supports saving state to file.
	string getState();
	//! set QGLViewer state from string (XML); QGLVIewer normally only supports loading state from file.
	void setState(string);
	//! Load display parameters (QGLViewer and OpenGLRenderer) from Scene::dispParams[n] and use them
	void useDisplayParameters(size_t n);
	//! Save display parameters (QGOViewer and OpenGLRenderer) to Scene::dispParams[n]
	void saveDisplayParameters(size_t n);
	//! Get radius & center of the part of scene that fits the current view
	std::pair<double, qglviewer::Vec> displayedSceneRadiusCenter();

	//! Adds our attributes to the QGLViewer state that can be saved
	QDomElement domElement(const QString& name, QDomDocument& document) const override;
	//! Adds our attributes to the QGLViewer state that can be restored
	void initFromDOMElement(const QDomElement& element) override;

	// if defined, snapshot will be saved to this file right after being drawn and the string will be reset.
	// this way the caller will be notified of the frame being saved successfully.
	string nextFrameSnapshotFilename;
#ifdef YADE_GL2PS
	// output stream for gl2ps; initialized as needed
	FILE* gl2psStream;
#endif

	boost::posix_time::ptime getLastUserEvent() const;


	DECLARE_LOGGER;

protected:
	void keyPressEvent(QKeyEvent* e) override;
	void postDraw() override;
	// overridden in the player that doesn't get time from system clock but from the db
	virtual string getRealTimeString() const;
	void           closeEvent(QCloseEvent* e) override;
	void           postSelection(const QPoint& point) override;
	void           endSelection(const QPoint& point) override;
	void           mouseDoubleClickEvent(QMouseEvent* e) override;
	void           wheelEvent(QWheelEvent* e) override;
	void           mouseMoveEvent(QMouseEvent* e) override;
	void           mousePressEvent(QMouseEvent* e) override;
	QString        helpString() const override;

	// Draws text, where each letter has a shifted background letter of opposite color.
	void drawReadableNum(const Real& n, const Vector3r& pos, unsigned precision = 4, const Vector3r& color = Vector3r(1, 1, 1));
	void drawReadableText(const std::string& txt, const Vector3r& pos, const Vector3r& color = Vector3r(1, 1, 1));
	void drawTextWithPixelShift(const std::string& txt, const Vector3r& pos, const Vector2i& shift, const Vector3r& color = Vector3r(1, 1, 1));
};

class YadeCamera : public qglviewer::Camera {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wsuggest-override"
	Q_OBJECT
#pragma GCC diagnostic pop
private:
	float cuttingDistance;

public:
	YadeCamera()
	        : cuttingDistance(0) {};
#if QGLVIEWER_VERSION >= 0x020603
	using QGLCompatDouble = qreal;
#else
	using QGLCompatDouble = float;
#endif
	QGLCompatDouble zNear() const override;

	virtual void setCuttingDistance(float s) { cuttingDistance = s; };
};

} // namespace yade
