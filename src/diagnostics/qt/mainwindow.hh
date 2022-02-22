//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/qt/mainwindow.hh
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Declarations for draco info main Qt window.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef diagnostics_qt_mainwindow_hh
#define diagnostics_qt_mainwindow_hh

#include "diWidget.hh"
#include <QMainWindow>

//namespace Ui {
//   class MainWindow;
//}

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = nullptr);
  ~MainWindow() override = default;

  // disable copy/move construction
  MainWindow(MainWindow const &rhs) = delete;
  MainWindow(MainWindow const &&rhs) = delete;

  // disable assignment and move-assignment
  MainWindow &operator=(MainWindow const &rhs) = delete;
  MainWindow &operator=(MainWindow const &&hs) = delete;

private slots:
  // None

private:
  // Forms
  //Ui::MainWindow *ui;

  // Widgets
  diWidget *diw;
};

#endif // diagnostics_qt_mainwindow_hh

//------------------------------------------------------------------------------------------------//
// end of diagnostics/qt/mainwindow.hh
//------------------------------------------------------------------------------------------------//
