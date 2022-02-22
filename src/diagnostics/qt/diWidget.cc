//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/qt/diWidget.cc
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Implementation for draco info widget.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "diWidget.hh"
#include "../draco_info.hh"

diWidget::diWidget(QWidget *parent)
    : QWidget(parent), layout(new QGridLayout(this)), pushbutton1(new QPushButton("&Ok")) {
  // Set the window title.
  setWindowTitle("draco_info-gui");

  // layout = new QGridLayout(this);

  // pushbutton1 = new QPushButton("&Ok");
  connect(pushbutton1, SIGNAL(clicked()), parent, SLOT(close()));

  // Get the release information from Lib_diagnostics
  QString msg(rtt_diagnostics::DracoInfo().briefReport().c_str());

  label1 = new QLabel(msg, this);
  label1->setWordWrap(false);

  layout->addWidget(label1, 0, 0);
  layout->addWidget(pushbutton1, 1, 0, Qt::AlignCenter);
  layout->columnStretch(0);

  setLayout(layout);
}

//------------------------------------------------------------------------------------------------//
// end of diagnostics/qt/diWidget.cc
//------------------------------------------------------------------------------------------------//
