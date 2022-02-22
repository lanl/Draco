//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/qt/mainwindow.cc
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Implementation for draco info main Qt window.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "mainwindow.hh"
#include "ui_mainwindow.h"
#include <QStatusBar>

//------------------------------------------------------------------------------------------------//
//! Constructor
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
      //ui(new Ui::MainWindow),
      diw(new diWidget(this)) {
  //ui->setupUi(this);
  setCentralWidget(diw);
  statusBar()->showMessage(tr("Status Bar"));
}

//------------------------------------------------------------------------------------------------//
// end of diagnostics/qt/mainwindow.cc
//------------------------------------------------------------------------------------------------//
