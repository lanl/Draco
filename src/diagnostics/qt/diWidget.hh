//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   diagnostics/qt/diWidget.hh
 * \author Kelly Thompson
 * \date   Monday, Aug 11, 2016, 17:05 pm
 * \brief  Declaration of draco info widget.
 * \note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef diwidget_hh
#define diwidget_hh

#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QWidget>

class diWidget : public QWidget {
  Q_OBJECT

public:
  explicit diWidget(QWidget *parent = nullptr);
  ~diWidget() override = default;

  // disable copy/move construction
  diWidget(diWidget const &rhs) = delete;
  diWidget(diWidget const &&rhs) = delete;

  // disable assignment and move-assignment
  diWidget &operator=(diWidget const &rhs) = delete;
  diWidget &operator=(diWidget const &&rhs) = delete;

private slots:
  /* void on_actionAbout_triggered(); */
  /* void on_actionDialog_with_Layout_triggered(); */
  /* void on_actionEdit_File_triggered(); */

private:
  QGridLayout *layout;
  QLabel *label1{nullptr};
  QPushButton *pushbutton1;
};

#endif // diwidget_hh

//------------------------------------------------------------------------------------------------//
// end of diagnostics/qt/diWidget.hh
//------------------------------------------------------------------------------------------------//
