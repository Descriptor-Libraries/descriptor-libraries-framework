
import { Route, Routes } from 'react-router-dom';
import MoleculeInfo from '../pages/Molecule';


export const AppRoutes = ({ pages }) => (
  <Routes>
    <Route path='/' element={pages['Home']}></Route>
    {Object.keys(pages).map((page, index) => (
      <Route key={index} path={`/${page.replace(" ", "_").toLowerCase()}`} element={pages[page]} />
    ))}
    <Route path="molecule/:molid" element={<MoleculeInfo />} />
  </Routes>
);